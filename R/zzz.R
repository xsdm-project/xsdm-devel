.onLoad <- function(libname, pkgname) {
  # Set future.globals.maxSize to your desired default value (e.g., 1 GB)
  options(future.globals.maxSize = 8.0 * 1024^3) # 1 GB in bytes
}
