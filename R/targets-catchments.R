#' Copies NHD files with `vpu` to `new_dir`.
#'
#' @param vpu VPU identifier
#' @param files List of files to scan
#' @param new_dir Directory to copy matching files into
#' @returns Path(s) of copied file(s)
#'
#' @keywords internal
rf.helper.copy_file <- function(vpu, files, new_dir) {
  f <- grep(vpu, files)
  tmp01 <- files[f]
  tmp02 <- file.path(new_dir, vpu, basename(files[f]))
  rf.utils.ensure_directory(dirname(tmp02))

  if (!file.exists(tmp02)) {
    cat("Copying", vpu, "to", paste0("'", new_dir, "'"))
    file.copy(grep(tmp01, files, value = TRUE), tmp02)
  }

  tmp02
}

#' @export
rf.targets.catchment_info <- function(vpu_topology, catchments_paths, flowlines_paths, dir_cleaned) {
  vpus <- unique(unlist(vpu_topology))

  cats <- data.frame(
    cat_path = catchments_paths,
    vpu = rf.utils.extract_vpu(catchments_paths)
  )

  fl <- data.frame(
    fl_path = flowlines_paths,
    vpu = rf.utils.extract_vpu(flowlines_paths)
  )

  dplyr::left_join(cats, fl, by = "vpu") |>
    dplyr::mutate(outfile = file.path(dir_cleaned, vpu, "catchments.fgb")) |>
    dplyr::filter(vpu %in% vpus) |>
    dplyr::relocate(vpu) |>
    dplyr::rowwise()
}

#' Clean Catchments Target
#' 
#' @description
#' Processes catchments according to steps 4.2.1 - 4.2.3.
#' 
#' @param cat_info `data.frame` containing path information for catchments and flowlines
#' @param simplify_keep Tolerance value for weighted Visvalingam simplification
#' @returns Path to cleaned catchment output (`cat_info$outfile`).
#' 
#' @export
rf.targets.clean_catchment <- function(cat_info, simplify_keep = 0.20) {

  # Transform catchments to Conus Albers
  catchment <-
    sf::read_sf(cat_info$cat_path) |>
    sf::st_transform(5070)
  names(catchment) <- tolower(names(catchment))

  # Transform flowlines to Conus Albers
  flowlines <-
    sf::read_sf(cat_info$fl_path) |>
    sf::st_transform(5070)
  names(flowlines) <- tolower(names(flowlines))

  old_tmpdir <- getOption("ms_tempdir")
  new_tmpdir <- file.path(dirname(cat_info$outfile), "tmp")

  rf.utils.ensure_directory(new_tmpdir)

  options(ms_tempdir = new_tmpdir)
  on.exit(unlink(new_tmpdir, recursive = TRUE, force = TRUE))

  # Clean fragments in catchments. This function handles
  # fragments and topology-preserving simplification,
  # as described in steps 4.2.1 - 4.2.3.
  out <- hydrofab::clean_geometry(
    catchments = catchment,
    flowlines = flowlines,
    ID = "featureid",
    fl_ID = "comid",
    crs = 5070,
    keep = simplify_keep,
    force = TRUE,
    sys = TRUE
  )

  if (!is.null(old_tmpdir)) {
    options(ms_tempdir = old_tmpdir)
  }

  # Output cleaned catchments
  rf.utils.ensure_directory(dirname(cat_info$outfile))
  sf::st_write(out, cat_info$outfile, layer = "catchments", quiet = TRUE, delete_dsn = TRUE)
  cat_info$outfile
}

#' Catchment VPU Rectification Target
#' 
#' @description
#' Processes VPU rectification of _cleaned_ catchments according to step 4.2.4.
#' 
#' @note
#' This process is highly data-dependent, so we must process it sequentially.
#' 
#' @param cat_cleaned_paths Paths to cleaned catchment files
#' @param vpu_topology `data.frame` containing topology definition for VPUs
#' @param dir_reference Output directory for reference features
#' 
#' @returns List of catchment output paths
#' 
#' @export
rf.targets.rectify_catchment_borders <- function(cat_cleaned_paths, vpu_topology, dir_reference) {
  already_processed <- list.files(dir_reference, pattern = "catchments\\.fgb", recursive = TRUE, full.names = TRUE)
  if (length(already_processed) == length(unique(unlist(vpu_topology)))) {
    return(already_processed)
  }

  tmpdir <- file.path(dir_reference, "tmp")
  rf.utils.ensure_directory(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE))

  # For each row of the VPU topology, we read in the catchments
  # corresponding to VPU A and VPU B. Then, we process them along
  # the boundary, and write the catchments back to their corresponding
  # output file.
  for (i in seq_len(nrow(vpu_topology))) {
    VPU1 <- vpu_topology$VPU1[i] # A
    VPU2 <- vpu_topology$VPU2[i] # B

    v_path_1 <-
      rf.helper.copy_file(VPU1, cat_cleaned_paths, dir_reference)
    v_path_2 <-
      rf.helper.copy_file(VPU2, cat_cleaned_paths, dir_reference)

    # Read in catchments for VPU A/B. We transform to 4326 to ensure
    # consistency with mapshaper. (Otherwise, it complains about planar SRS).
    vpu_div_1 <-
      sf::read_sf(v_path_1, "catchments") |>
      dplyr::mutate(vpuid = VPU1) |>
      sf::st_transform(4326)

    vpu_div_2 <-
      sf::read_sf(v_path_2, "catchments") |>
      dplyr::mutate(vpuid = VPU2) |>
      sf::st_transform(4326)

    # Join VPU A and VPU B into a single file for mapshaper to clean
    tmpfile <- tempfile(pattern = "rectify_borders", tmpdir = tmpdir, fileext = ".geojson")
    dplyr::bind_rows(vpu_div_1, vpu_div_2) |>
      yyjsonr::write_geojson_file(tmpfile)

    system2("mapshaper", args = c(tmpfile, "-clean", "-o", "force", tmpfile))

    # Read in cleaned catchment files and ensure area is consistent
    new <-
      yyjsonr::read_geojson_file(tmpfile) |>
      sf::st_as_sf() |>
      sf::st_set_crs(4326) |>
      sf::st_transform(5070) |>
      dplyr::select(featureid, vpuid, geometry) |>
      dplyr::mutate(areasqkm = hydrofab::add_areasqkm(geometry))

    unlink(tmpfile)

    # separate and write back VPU A
    dplyr::filter(new, vpuid == VPU1) |>
      sf::st_write(v_path_1, "catchments", delete_dsn = TRUE, quiet = TRUE)

    # separate and write back VPU B
    dplyr::filter(new, vpuid == VPU2) |>
        sf::st_write(v_path_2, "catchments", delete_dsn = TRUE, quiet = TRUE)
  }

  list.files(dir_reference, pattern = "catchments\\.fgb", recursive = TRUE, full.names = TRUE)
}

#' Catchments Reference Output Target
#' 
#' @description
#' Processes remaining catchment steps according to step 4.2.5.
#' 
#' @param rf_cat_path Path to VPU-rectified catchment
#' @returns Path to outputted reference catchment
#' 
#' @section FIXME: Snap to underlying grid with size of .0009?
#' 
#' @export
rf.targets.reference_catchments <- function(rf_cat_path) {

  cats <- sf::read_sf(rf_cat_path)

  # Handle catchments fully within other catchments
  imap <- sf::st_within(cats)

  df <- data.frame(
    within = rep(cats$featureid, times = lengths(imap)),
    featureid = cats$featureid[unlist(imap)]
  ) |>
    dplyr::filter(featureid != within)

  rm(imap)

  cats2 <- dplyr::filter(cats, !featureid %in% unlist(df))

  if (nrow(df) > 0) {
    # For each consuming catchment, combined the contained catchments and erase.
    d <- lapply(unique(df$featureid), FUN = function(id) {
      dplyr::filter(cats, featureid %in% dplyr::filter(df, featureid == id)$within) |>
        sf::st_combine() |>
        sf::st_make_valid() |>
        sf::st_difference(x = dplyr::filter(cats, featureid == id))
    })

    d <- sf::st_as_sf(dplyr::bind_rows(d))

    cats2 <- dplyr::bind_rows(cats2, d)

    rm(d)
  }

  # Combine differenced catchments, recalc area, and output to disk
  dplyr::filter(
    cats,
    featureid %in% dplyr::filter(df, !within %in% cats2$featureid)$within
) |>
    dplyr::bind_rows(cats2) |>
    sf::st_cast("POLYGON") |>
    dplyr::mutate(areasqkm = hydrofab::add_areasqkm(geometry)) |>
    sf::st_write(rf_cat_path, "catchments", delete_dsn = TRUE, quiet = TRUE)

  rf_cat_path
}
