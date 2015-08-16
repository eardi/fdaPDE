.noGenerics <- TRUE

.onLoad <- function(libname, pkgname) 
{
  library.dynam("FEMr", pkgname, libname)
}

.onUnload <- function(libpath)
  library.dynam.unload("FEMr", libpath)
