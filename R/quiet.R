#' Suppress Console Messages and Output
#'
#' This function suppresses messages and output from expressions, optionally allowing message output to display.
#'
#' @param ... Expressions to evaluate quietly.
#' @param messages Logical; if `TRUE`, displays messages. If `FALSE` (default), suppresses messages.
#' @param cat Logical; if `TRUE`, displays console output. If `FALSE` (default), suppresses console output by redirecting to a temporary file.
#'
#' @return The result of the evaluated expressions in `...`, with messages or output suppressed based on `messages` and `cat`.
#' @noRd
#'
#' @examples
#' # Example: Run a command quietly without messages or output
#' quiet({message("This is a test message")}, messages = FALSE, cat = FALSE)
#'
#' # Example: Allow messages but suppress other output
#' quiet({message("This will show"); print("This will not")}, messages = TRUE, cat = FALSE)
quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    tmpf <- tempfile()
    sink(tmpf)
    on.exit({sink(); file.remove(tmpf)})
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}
