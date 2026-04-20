xsdmStartupMessage <- function()
{
  msg <- c(paste0(
    "              _
__  _____  __| |_ __ ___
\\ \\/ / __|/ _` | '_ ` _ \\
 >  <\\__ \\ (_| | | | | | |
/_/\\_\\___/\\__,_|_| |_| |_|    version ",
    utils::packageVersion("xsdm")),
    "\nType 'citation(\"xsdm\")' for citing this R package in publications.\n",
    "\nThis package depends on cmdstanr that is not in CRAN\n"
  )
  
  return(msg)
}

.onLoad <- function(libname, pkgname) {
  # Set future.globals.maxSize to your desired default value (e.g., 1 GB)
  options(future.globals.maxSize = 8.0 * 1024^3) # 1 GB in bytes
}

.onAttach <- function(lib, pkg)
{
  # unlock .xsdm variable allowing its modification
  unlockBinding("xsdm", asNamespace("xsdm"))
  # startup message
  msg <- xsdmStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'xsdm' version", utils::packageVersion("xsdm"))
  packageStartupMessage(msg)
  invisible()
  
}
.onUnload <- function(libpath) {
  library.dynam.unload("xsdm", libpath)
}