#' A simple wrapper for the SCOUP executable
#'
#' @param method which executable to run, can be one of \code{"scoup"}, \code{"scoup_resume"}, \code{"sp"}, or \code{"correlation"}.
#' @param args a string containing the command that would normally be
#'   entered for running the SCOUP executable.
#' @param verbose whether or not to print the SCOUP output
#'
#' @export
run_SCOUP <- function(method, args_string, verbose = F) {
  args <- strsplit(args_string, " ")[[1]]
  method_int <- c("scoup"=1, "scoup_resume"=2, "correlation"=3, "sp"=4)[method]
  main_wrap(method_int, args, verbose)
}
