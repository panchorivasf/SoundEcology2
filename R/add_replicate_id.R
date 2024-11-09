#' Add replicate id for a dataset with severall treatments 
#'
#' @param data a data frame 
#' @param treatment the name of the column containing the treatment identifiers
#'
#' @return a data frame with treatment_id column
#' @export
#' 
#' @importFrom dplyr mutate mutate
#'
#' @examples add_replicate_id(my_dataset, treatment = 'zone')
add_replicate_id <- function(data, treatment = "treatment") {
  data |>
    group_by(treatment) |> 
    mutate(replicate = dense_rank(sensor_id)) |> 
    ungroup()
}