.noGenerics <- TRUE

.onLoad <- function(libname, pkgname) 
{
  library.dynam("rFEM", pkgname, libname)
}

.onUnload <- function(libpath)
  library.dynam.unload("rFEM", libpath)
