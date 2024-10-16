#' Burnline Events Processing Target
#' 
#' @description
#' Processes burn line events according to step 4.3.1
#' 
#' @param nhd_gdb Path to NHDPlus Seamless GeoDatabase
#' @param rf_vaa NHD Value-added attributes table
#' @param dir_data Data output directory
#' 
#' @returns Path to outputted burn line events file for CONUS
#' 
#' @export
rf.targets.burnline_events <- function(nhd_gdb, rf_vaa, dir_data) {
  ble_outfile <- file.path(dir_data, "burnline_events.fgb")

  # Identify all flowlines that are a headwater (startflag == 1) or
  # a minor path from divergence (divergence == 2).
  nhdflag <-
    rf_vaa |>
    dplyr::filter(startflag == 1 | divergence == 2) |>
    dplyr::select(COMID = comid, vpuid)

  sf::read_sf(nhd_gdb, query = "SELECT COMID, Shape AS geometry FROM BurnLineEvent") |>
    # Drop Z and M dimensions
    sf::st_zm() |>
    dplyr::right_join(nhdflag, by = "COMID") |>
    sf::st_as_sf() |>
    # Remove empty geometries
    dplyr::filter(!sf::st_is_empty(geometry)) |>
    dplyr::select(COMID, vpuid, geometry) |>
    sf::st_cast("LINESTRING") |>
    sf::write_sf(ble_outfile)

  rm(nhdflag)

  ble_outfile
}
