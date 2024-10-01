#' Print a list vertically
#'
#' @param list a list object.
#'
#' @return a printed list.
#' @export
#'
#' @examples print_list(list)
print_list <- function(list){
  for (i in seq_along(list)) {
    cat(i, list[[i]], "\n")
  }
}