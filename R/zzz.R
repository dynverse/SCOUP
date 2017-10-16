.onLoad <- function(libname, pkgname) {
  if (!file.exists(paste0(find.package("SCOUP"), "/code/scoup"))) {
    reinstall()
  }
}

reinstall <- function() {
  path <- find.package("SCOUP")
  cmd <- paste0("cd '", path, "/code'; make all")
  system(cmd)
}
