#' Clean Flowlines Target
#' 
#' @description
#' Processes flowlines according to step 4.3.2
#' 
#' @param flowlines_path Path to NHD flowlines file
#' @param rf_cat_paths Path to _reference_ catchment outputs
#' @param rf_ble_path Path to CONUS burn line events file
#' @param rf_vaa NHDPlus-VAA table
#' @param rf_enhd E2NHDPlus River Attributes table
#' @param dir_cleaned Cleaned output directory
#' @returns Path to cleaned flowlines
#' 
#' @export
rf.targets.clean_flowlines <- function(flowlines_path, rf_cat_paths, rf_ble_path, rf_vaa, rf_enhd, dir_cleaned) {
  nhd <-
    sf::read_sf(flowlines_path) |>
    sf::st_zm() |>
    sf::st_transform(5070) |>
    nhdplusTools::align_nhdplus_names() |>
    dplyr::left_join(rf_vaa, by = c("COMID" = "comid")) |>
    dplyr::select(
      COMID,
      fromnode,
      tonode,
      startflag,
      streamcalc,
      divergence,
      dnminorhyd
    ) |>
    dplyr::left_join(rf_enhd, by = c("COMID" = "comid")) |>
    nhdplusTools::align_nhdplus_names() |>
    dplyr::mutate(LENGTHKM = hydrofab::add_lengthkm(geometry))

  rm(rf_vaa)

  cats <-
    rf_cat_paths[grep(pattern = rf.utils.extract_vpu(flowlines_path), x = basename(dirname(rf_cat_paths)))] |>
    sf::read_sf() |>
    sf::st_transform(5070)

  rm(rf_cat_paths)

  # Here we begin checking if a flowline's starting node lies in a catchemnt
  # other than the one it is associated with. If it is, we replace its geometry
  # with the BLE geometry (reducing the effective length by 150m).

  ble_comid <- sf::read_sf(rf_ble_path, query = "SELECT COMID FROM burnline_events")$COMID
  
  # Get the flowlines with a matching burn line event
  ble_options <-
    dplyr::filter(nhd, COMID %in% ble_comid) |>
    sf::st_transform(5070)

  rm(ble_comid)
  
  ## Convert matching flowlines' geometry to their starting node
  sf::st_geometry(ble_options) <-
    nhdplusTools::get_node(ble_options, "start") |>
    sf::st_geometry()
  
  # Check starting node
  imap <- sf::st_intersects(ble_options, cats)

  df <- data.frame(
    start = rep(ble_options$COMID, times = lengths(imap)),
    cat = cats$featureid[unlist(imap)]
  ) |>
    dplyr::filter(start != cat)

  rm(cats)
  rm(ble_options)
  rm(imap)

  # Get corresponding burn line geometries
  ble <- sf::read_sf(rf_ble_path, query = paste0(
    "SELECT COMID, '_ogr_geometry_' AS geometry FROM burnline_events WHERE vpuid = '", rf.utils.extract_vpu(flowlines_path), "'"
  )) |>
    dplyr::filter(COMID %in% df$start) |>
    dplyr::distinct(COMID, .keep_all = TRUE) |>
    sf::st_as_sf() |>
    sf::st_transform(5070)

  # Indices for which burnline matches what flowline
  .matches <- match(ble$COMID, nhd$COMID)

  rm(df)

  # This is where we replace the geometry
  sf::st_geometry(nhd)[.matches] <- sf::st_geometry(ble)
  nhd$ble <- FALSE
  nhd$ble[.matches] <- TRUE

  rm(ble)

  custom_net <-
    nhd |>
    sf::st_drop_geometry() |>
    dplyr::select(COMID, FromNode, ToNode, Divergence) |>
    nhdplusTools::get_tocomid(remove_coastal = FALSE) |>
    dplyr::select(comid, override_tocomid = tocomid)

  outfile <- file.path(dir_cleaned, rf.utils.extract_vpu(flowlines_path), "flowlines.fgb")

  nhd |>
    dplyr::left_join(custom_net, by = c("COMID" = "comid")) |>
    dplyr::mutate(override_tocomid = ifelse(toCOMID == 0, override_tocomid, toCOMID)) |>
    sf::st_write(outfile, "flowlines", delete_dsn = TRUE, quiet = TRUE)

  outfile
}

#' Reference Flowlines Target
#' 
#' @description
#' Processes flowlines according to step 4.3.3
#' 
#' @param nhd Path to _cleaned_ flowlines
#' @param vpu Current VPU for `nhd`.
#' @param rf_enhd_comid COMIDs with E2NHDPlus river attributes
#' @param dir_reference Reference output directory
#' @returns Path to reference flowlines
#' 
#' @export
rf.targets.reference_flowlines <- function(nhd, vpu, rf_enhd_comid, dir_reference) {
  nhd <- sf::read_sf(nhd)

  check <- !nhd$COMID %in% nhd$override_tocomid & !(
    nhd$override_tocomid == 0 |
      is.na(nhd$override_tocomid) |
      !nhd$override_tocomid %in% nhd$COMID
  )

  check_direction <- dplyr::filter(nhd, check)

  if (!all(
    check_direction$override_tocomid[
      check_direction$override_tocomid != 0
    ] %in% nhd$COMID
  )) {
    targets::tar_throw_validate(
      "Not all of non-zero override_tocomid appear in nhd$COMID"
    )
  }

  .matches <- match(check_direction$override_tocomid, nhd$COMID)
  .matches <- nhd[.matches, ]

  fn_list <- Map(
    list,
    split(
      check_direction,
      seq_len(nrow(check_direction))
    ),
    split(
      .matches,
      seq_len(nrow(.matches))
    )
  )

  new_geom <- lapply(fn_list, FUN = function(x) {
    nhdplusTools::fix_flowdir(
      x[[1]]$COMID,
      fn_list = list(
        flowline  = x[[1]],
        network   = x[[2]],
        check_end = "end"
      )
    )
  })

  error_index <- sapply(new_geom, inherits, what = "try-error")

  if (length(error_index) == 0) {
    # Prevent class(error_index) == list, which causes an error
    errors <- head(nhd, 0)
  } else {
    errors <- dplyr::filter(
      nhd,
      COMID %in% check_direction$COMID[error_index]
    )
  }

  check[which(nhd$COMID %in% errors$COMID)] <- FALSE

  if (!all(sapply(sf::st_geometry(errors), sf::st_is_empty))) {
    targets::tar_throw_validate(
      "Errors other than empty geometry from fix flowdir"
    )
  }

  ng <- do.call(c, new_geom[!error_index])

  ng_direction <- ng |>
    (\(geom) lwgeom::st_endpoint(geom) - lwgeom::st_startpoint(geom))() |>
    sf::st_coordinates()
  colnames(ng_direction) <- c("X1", "Y1")

  nhd_direction <- sf::st_geometry(nhd)[check] |>
    (\(geom) lwgeom::st_endpoint(geom) - lwgeom::st_startpoint(geom))() |>
    sf::st_coordinates()
  colnames(nhd_direction) <- c("X2", "Y2")

  switched <-
    cbind(ng_direction, -nhd_direction) |>
    rowSums() |>
    abs() |>
    (\(x) x >= 1e-5)()

  nhd$reversed <- FALSE
  nhd$reversed[check] <- switched
  sf::st_geometry(nhd)[check] <- ng

  outfile <- file.path(dir_reference, vpu, "flowlines.fgb")
  rf.utils.ensure_directory(dirname(outfile))

  nhd |>
    dplyr::filter(COMID %in% rf_enhd_comid) |>
    dplyr::select(-override_tocomid, -LENGTHKM) |>
    dplyr::mutate(lengthkm = hydrofab::add_lengthkm(geometry)) |>
    sf::st_cast("LINESTRING") |>
    sf::st_write(outfile, "flowlines", delete_dsn = TRUE, quiet = TRUE)

  outfile
}
