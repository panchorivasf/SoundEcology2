#' Extract and Display Function Arguments and Default Values
#'
#' This function extracts the arguments of a given function along with their
#' default values and displays them in a formatted vertical layout on the console.
#' It also returns a named list containing the argument names and their defaults.
#'
#' @param func A function object whose arguments and defaults should be extracted.
#'   Can be any R function including built-in functions, user-defined functions,
#'   or functions from loaded packages.
#'
#' @return A named list where names are the argument names and values are the
#'   default values. Arguments without defaults are marked as "No default".
#'   The list is returned invisibly to avoid cluttering the console output.
#'
#' @details
#' The function handles different types of default values:
#' \itemize{
#'   \item Arguments with no defaults (marked as "No default")
#'   \item Arguments with literal values (numbers, strings, logical values)
#'   \item Arguments with expressions or function calls as defaults
#' }
#'
#' The console output is formatted with argument names left-aligned and
#' default values clearly displayed. For functions with no arguments,
#' an appropriate message is shown.
#'
#' @examples
#' # Extract defaults from built-in functions
#' defaults(lm)
#' defaults(plot)
#'
#' # Extract defaults from a custom function
#' my_func <- function(x, y = 10, z = "test", data = NULL) {
#'   return(x + y)
#' }
#' result <- defaults(my_func)
#'
#' # The returned list can be used programmatically
#' args_info <- defaults(mean)
#' print(args_info)
#'
#' @seealso
#' \code{\link{formals}} for extracting function arguments,
#' \code{\link{args}} for a simpler argument display
#'
#' @export
#'
#' @author Your Name
#'