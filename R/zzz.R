.onLoad <- function(libname, pkgname) {
  options(future.globals.maxSize = 8.0 * 1024^3)
  assign("xsdm", new.env(), envir = parent.env(environment()))
}

xsdmStartupMessage <- function()
{
  msg <- c(paste0(
    "              _
__  _____  __| |_ __ ___
\\ \\/ / __|/ _` | '_ ` _ \\
 >  <\\__ \\ (_| | | | | | |
/_/\\_\\___/\\__,_|_| |_| |_|    version ",
    utils::packageVersion("xsdm")),
    "\nType 'citation(\"xsdm\")' for citing this R package in publications.\n"
  )
  
  return(msg)
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
