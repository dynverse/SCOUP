#' #' A simple wrapper for the GEDEVO executables
#' #'
#' #' @param args a string containing the command that would normally be
#' #'   entered for running the GEDEVO executable.
#' #' @param verbose whether or not to print the GEDEVO output
#' #'
#' #' @examples
#' #' \dontrun{
#' #'   cmd <- paste0(
#' #'     "./gedevo --resume test.state --create --groups DM hprd --ntw DM.ntw "
#' #'     "--ntw hprd.ntw --grsig hprd.sigs hprd --grsig DM.sigs DM --undirected "
#' #'     "--pop 1000 --maxsame 500 --no-prematch")
#' #'   run_GEDEVO(cmd)
#' #' }
#' #'
#' #' @export
#' run_GEDEVO <- function(args_string, verbose = F) {
#'   args <- strsplit(args_string, " ")[[1]]
#'   poor_mans_wrap(args, verbose)
#' }
