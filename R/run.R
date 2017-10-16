execute <- function(method, args_string, verbose = FALSE) {
  path <- find.package("SCOUP")
  cmd <- paste0("'", path, "/code/", method, "' ", args_string)
  system(cmd, ignore.stdout = !verbose, ignore.stderr = !verbose)
}

#' Run SCOUP
#'
#' @param expr The expression data
#' @param start_ix The indices or names of the starting population
#' @param ndim ndim
#' @param nbranch nbranch
#' @param max_ite1 max_ite1
#' @param max_ite2 max_ite2
#' @param alpha_min alpha_min
#' @param alpha_max alpha_max
#' @param t_min t_min
#' @param t_max t_max
#' @param sigma_squared_min sigma_squared_min
#' @param thresh thresh
#'
#' @importFrom stats var
#' @importFrom utils write.table read.table
#' @importFrom glue glue
#'
#' @export
run_SCOUP <- function(expr,
                      start_ix,
                      ndim = 3,
                      nbranch = 1,
                      max_ite1 = 1000,
                      max_ite2 = 10000,
                      alpha_min = .1,
                      alpha_max = 100,
                      t_min = .001,
                      t_max = 2,
                      sigma_squared_min = .1,
                      thresh = .01,
                      verbose = FALSE) {
  # create distribution on starting population
  vars <- apply(expr[start_ix,], 2, stats::var)
  means <- apply(expr[start_ix,], 2, mean)
  distr_df <- data.frame(i = seq_along(vars) - 1, means, vars)

  # create temporary folder
  tmp_dir <- tempfile(pattern = "scoup")
  dir.create(tmp_dir)

  tryCatch({
    # write data to files
    utils::write.table(t(expr), file = paste0(tmp_dir, "/data"), sep = "\t", row.names = FALSE, col.names = FALSE)
    utils::write.table(distr_df, file = paste0(tmp_dir, "/init"), sep = "\t", row.names = FALSE, col.names = FALSE)

    # execute sp
    execute("sp", glue::glue(
      "{tmp_dir}/data",
      "{tmp_dir}/init",
      "{tmp_dir}/time_sp",
      "{tmp_dir}/dimred",
      "{ncol(expr)}",
      "{nrow(expr)}",
      "{ndim}",
      .sep = " "
    ), verbose = verbose)

    # execute scoup
    execute("scoup", glue::glue(
      "{tmp_dir}/data",
      "{tmp_dir}/init",
      "{tmp_dir}/time_sp",
      "{tmp_dir}/gpara",
      "{tmp_dir}/cpara",
      "{tmp_dir}/ll",
      "{ncol(expr)}",
      "{nrow(expr)}",
      "-k {nbranch}",
      "-m {max_ite1}",
      "-M {max_ite2}",
      "-a {alpha_min}",
      "-A {alpha_max}",
      "-t {t_min}",
      "-T {t_max}",
      "-s {sigma_squared_min}",
      "-e {thresh}",
      .sep = " "
    ), verbose = verbose)

    # read dimred
    dimred <- utils::read.table(
      paste0(tmp_dir, "/dimred"),
      col.names = c("i", paste0("Comp", seq_len(ndim)))
    )

    # last line is root node
    root <- dimred[nrow(dimred),-1,drop=F]
    dimred <- dimred[-nrow(dimred),]
    rownames(dimred) <- rownames(expr)

    # read cell params
    cpara <- utils::read.table(
      paste0(tmp_dir, "/cpara"),
      col.names = c("time", paste0("M", seq_len(nbranch)))
    )
    rownames(cpara) <- rownames(expr)

    # read gene params
    pi <- utils::read.table(paste0(tmp_dir, "/gpara"), nrows = 1)[1,]
    gpara <- utils::read.table(
      paste0(tmp_dir, "/gpara"),
      skip = 1,
      col.names = c("alpha", "sigma", paste0("theta_", seq_len(nbranch)))
    )
    rownames(gpara) <- colnames(expr)

    # loglik
    ll <- utils::read.table(paste0(tmp_dir, "/ll"))[[1]]

    # return output
    list(
      root = root,
      dimred = dimred,
      cpara = cpara,
      pi = pi,
      gpara = gpara,
      ll = ll
    )
  }, finally = {
    # unlink(tmp_dir, recursive = TRUE, force = TRUE)
  })
}
