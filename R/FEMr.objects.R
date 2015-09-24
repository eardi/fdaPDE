#' Define a Piecewise Polinomial basis
#' 
#' @param mesh An object of class \code{TRIMESH2D} representing a triangulation. See \link{create.MESH.2D}.
#' @param order Order of the polynomial restricted to each element of the triangulaton, 
#' which may be either 1 or 2. When ORDER = 1 the basis system is piecewise linear. When 
#' ORDER = 2 the basis system is piecewise quadratic. This parameter must be less or equal respect to the order 
#' specified in the mesh creation through \link{create.MESH.2D}.
#' @param CPP_CODE if \code{TRUE} avoids the computation of some additional elements, 
#' not necessary if the functions depending on the created basis are called with the same flag 
#' \code{CPP_CODE=TRUE}
#' @return An object of class \code{FEM}. This contains the \code{mesh} along with the following variables:
#' \item{\code{order}}{Order of elements.} 
#' \item{\code{nbasis}}{The number of basis.} 
#' \item{\code{J}}{The area of each triangle of the basis.} 
#' \item{\code{transf}}{A matrix such that \code{transf[i,,]} is the 2-by-2 tranformation matrix that transforms the nodes of the reference triangle to the nodes of the i-th triangle.}
#' \item{\code{metric}}{A matrix such that \code{metric[i,,]} is the 2-by-2 matrix \code{transf[i,,]^{-1}*transf[i,,]^{-T}}. This matrix is usuful for the computation
#' of the integrals over the elements of the mesh. This is necessary for the computation of the discretized problem.}
#' @description Sets up a Finite Element basis for the analysis of functional bivariate data. This
#' defines the discetized functional space where the infinite-dimensional formulations are solved.
#' The basis' function are globally continuos, and linear or quadratic polynomial if restricted 
#' to each element of the triangulation. See e.g. Sangalli, Ramsay, Ramsay (2013).
#' @references Sangalli, L.M., Ramsay, J.O. & Ramsay, T.O., 2013. Spatial spline regression models. Journal of the Royal Statistical Society. Series B: Statistical Methodology, 75(4), pp.681-703.
#' @usage create.FEM.basis(mesh, order, CPP_CODE = FALSE)
#' @examples 
#' ## Creates an object TRIMESH2D with a concavity and second order nodes
#' mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
#' segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)), order=1)
#' ## Plot it
#' plot(mesh)                   
#' ## Creates the basis object
#' basisobj = create.FEM.basis(mesh, order = 1)

create.FEM.basis = function(mesh, order, CPP_CODE = FALSE)
{
  type <- NULL
  if (class(mesh) == "TRIMESH2D")
  {
    type <- "FEM"
  }
  
  #  The number of basis functions corresponds to the number of vertices
  #  for order = 1, and to vertices plus edge midpoints for order = 2
  
  nbasis = dim(mesh$nodes)[[1]]
  
  eleProp = NULL
  if(CPP_CODE == FALSE)
  {
    eleProp = R_elementProperties(mesh)
  }
  
  basisobj = list(type = type, mesh = mesh, order = as.integer(order), nbasis = nbasis, J=eleProp$J, transf = eleProp$transf, metric = eleProp$metric)
  
  basisobj
}


#' Define a Functional Data through a Finite Element basis
#' 
#' @param coefmat A vector or a matrix containg the coefficients of a Finite Element representation. The number of rows (or the vector's length) corresponds to the number of functions defining the basis, i.e.
#' the number of elements in the mesh. The number of columns corresponds to the number of functional replicates. 
#' @param basisobj An object defining the Finite Element basis. It can be created with \link{create.FEM.basis}.
#' @description This is the constructur function for objects of the class FEM. Each function of FEMr that needs to create a 
#' new Functional Object defined respect to the Finite Element (FE) basis, must call this function. Usually users do not need to call this function directly.
#' @usage FEM(coefmat,basisobj)
#' @return An object of the class \code{FEM}. This contains a list with components \code{coefmat} and \code{basis}.
#' @examples 
#' library(FEMr)
#' data("mesh.2D.rectangular")
#' plot(mesh.2D.rectangular)
#' ## FEM object with 1st order Finite Element basis
#' basisobj = create.FEM.basis(mesh.2D.rectangular, 1)
#' coeff <- sin(mesh.2D.rectangular$nodes[,1])*cos(mesh.2D.rectangular$nodes[,2])
#' FEM_object<- FEM(coeff, basisobj)
#' plot(FEM_object)

FEM<-function(coefmat,basisobj)
{
  if(is.vector(coefmat))
  {
    coefmat = as.matrix(coefmat)
  }
  fclass = NULL
  fclass = list(coefmat=coefmat, basisobj=basisobj)
  class(fclass)<-"FEM"
  return(fclass)
}

#' Plot a Functional Object of the class FEM
#' 
#' @param x A functional object of the class FEM.
#' @param num_refinements An natural number specifying how many bisections should by applied to each triangular element for
#' plotting purpose. This functionality is useful where a discretization with 2nd order Finite Element has been applied.
#' @param ... Arguments representing the graphical parameters to be passed to methods, see \link[rgl]{plot3d}.
#' @description Functional objects created with the function \code{FEM} or returned by \code{smooth.FEM.basis}, \code{smooth.FEM.PDE.basis} or
#' \code{smooth.FEM.PDE.SV.basis} can be plotted in 3D.
#' @usage \method{plot}{FEM}(x, num_refinements, ...)  
#' @examples 
#' data(MeuseData)
#' data(MeuseBorder)
#' order=1
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = order)
#' data = log(MeuseData[,7])
#' basisobj = create.FEM.basis(mesh, order)
#' lambda = 10^3.5
#' ZincMeuse = smooth.FEM.basis(observations = data, basisobj = basisobj, lambda = lambda)
#' ## plot of the FEM object representing the fitted function
#' plot(ZincMeuse$fit.FEM)
#' ## plot the FEM object representing the misfit
#' plot(ZincMeuse$PDEmisfit.FEM)

plot.FEM = function(x, num_refinements = NULL, ...)  
{
  if(x$basisobj$order == 1)
  {
    R_plot.ORD1.FEM(x, ...)
  }else{
    R_plot.ORDN.FEM(x, num_refinements, ...)
  }
}

#' Draw a Image Plot of a Functional Object of the class FEM
#' 
#' @param x A functional object of the class FEM.
#' @param num_refinements An natural number specifying how many bisections should by applied to each triangular element for
#' plotting purpose. This functionality is useful where a discretization with 2nd order Finite Element has been applied.
#' @param ... Arguments representing the graphical parameters to be passed to methods, see \link[rgl]{plot3d}.
#' @description Functional objects created with the function \code{FEM} or returned by \code{smooth.FEM.basis}, \code{smooth.FEM.PDE.basis} or
#' \code{smooth.FEM.PDE.SV.basis} can be visualized trough an image plot.
#' @usage \method{image}{FEM}(x, num_refinements, ...)  
#' @examples 
#' data(MeuseData)
#' data(MeuseBorder)
#' order=1
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = order)
#' data = log(MeuseData[,7])
#' basisobj = create.FEM.basis(mesh, order)
#' lambda = 10^3.5
#' ZincMeuse = smooth.FEM.basis(observations = data, basisobj = basisobj, lambda = lambda)
#' ## plot of the FEM object representing the fitted function
#' image(ZincMeuse$fit.FEM)
#' ## plot the FEM object representing the misfit
#' image(ZincMeuse$PDEmisfit.FEM)
image.FEM = function(x, num_refinements = NULL, ...)  
{
  if(x$basisobj$order == 1)
  {
    R_image.ORD1.FEM(x, ...)
  }else{
    R_image.ORDN.FEM(x, num_refinements, ...)
  }
}
