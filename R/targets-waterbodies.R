#' Clean Waterbodies Target
#' 
#' @description
#' Processes waterbodies according to steps 4.1.1 - 4.1.6.
#' 
#' @param waterbodies_path Path to unmodified NHDWaterbody file
#' @param dir_cleaned Directory to output cleaned waterbody geometry
#' @returns Path to cleaned waterbodies
#'
#' @export
rf.targets.clean_waterbodies <- function(waterbodies_path, dir_cleaned) {
  wb <-
    sf::read_sf(waterbodies_path) |>
    dplyr::rename_with(~ toupper(.x), .cols = c(dplyr::everything(), -geometry)) |>
    # Only keep features with FTYPE == LakePond or Reservoir
    dplyr::filter(FTYPE %in% c("LakePond", "Reservoir")) |>
    # Drop Z and M dimensions
    sf::st_zm() |>
    sf::st_make_valid()

  wb_clean <-
    # Dissolve features on COMID
    hydrofab::union_polygons(wb, "COMID") |>
    # Project to CONUS Albers Equal Area
    sf::st_transform(5070) |>
    # Remove Holes
    sfheaders::sf_remove_holes() |>
    # Ensure shapes are valid
    sf::st_make_value() |>
    dplyr::left_join(sf::st_drop_geometry(wb), by = "COMID") |>
    dplyr::select(GNIS_ID, GNIS_NAME, COMID, FTYPE, geometry)

  wb_outfile <- file.path(dir_cleaned, rf.utils.extract_vpu(waterbodies_path), "waterbodies.fgb")
  sf::st_write(wb_clean, wb_outfile, "waterbodies", delete_dsn = TRUE, quiet = TRUE)
  wb_outfile
}

#' Process Reference Waterbodies Target
#' 
#' @description
#' Processes _cleaned_ waterbodies according to steps 4.1.7 - 4.1.9.
#' 
#' @section FIXME: On/off network waterbodies based on NHD flowlines? 
#' 
#' @param wb_clean Path to cleaned waterbodies file
#' @param wb_vpu The corresponding VPU of `wb_clean`
#' @param wb_reference_dir Directory to output reference waterbodies
#' @returns Path to reference waterbodies
#'
#' @export
rf.targets.reference_waterbodies <- function(wb_clean, wb_vpu, wb_reference_dir) {
  wb_clean <-
    sf::read_sf(wb_clean)

  wb_filter <-
    wb_clean |>
      dplyr::filter(is.na(GNIS_ID)) |>
      dplyr::mutate(member_comid = as.character(COMID))
  
  # Implies that some GNIS_ID are not NA
  if (nrow(wb_filter) != nrow(wb_clean)) {
    wb_filter_2 <-
      wb_clean |>
      dplyr::filter(!is.na(GNIS_ID)) |>
      dplyr::mutate(area_sqkm = hydrofab::add_areasqkm(geometry))

    # Waterbodies with a GNIS_ID are dissolved on the GNIS_ID,
    # and added back to the set of waterbodies without a GNIS_ID.
    wb_filter <-
      wb_filter_2 |>
      hydrofab::union_polygons("GNIS_ID") |>
      dplyr::left_join(sf::st_drop_geometry(wb_filter_2), by = "GNIS_ID") |>
      dplyr::group_by(GNIS_ID) |>
      dplyr::mutate(member_comid = paste(COMID, collapse = ",")) |>
      dplyr::slice_max(area_sqkm) |>
      dplyr::ungroup() |>
      dplyr::bind_rows(wb_filter)

    rm(wb_filter_2)
  }

  wb_outfile <- file.path(wb_reference_dir, wb_vpu, "waterbodies.fgb")
  rf.utils.ensure_directory(dirname(wb_outfile))

  # Finalize and output reference waterbodies
  wb_filter |>
    dplyr::select(GNIS_ID, GNIS_NAME, COMID, FTYPE, member_comid, geometry) |>
    dplyr::mutate(area_sqkm = hydrofab::add_areasqkm(geometry)) |>
    sf::st_make_valid() |>
    sf::st_cast("POLYGON") |>
    sf::st_write(wb_outfile, "waterbodies", delete_dsn = TRUE, quiet = TRUE)

  wb_outfile
}
