#' Evaluate a FEM object on a set of point locations
#' 
#' @param locations A #locations-by-2 matrix where each row specifies the x and y coordinate of the corresponding location.
#' @param FEM the Functional Object of class \code{FEM} to be evaluated
#' @param CPP_CODE Boolean. If \code{TRUE} the computation relies on the C++ implementation of a Visibility Walk Algorithm (Devillers et al. 2001). This usually ensures a much faster computation.
#' @return 
#' A matrix of numeric evaluations of the \code{FEM} object. Each row indicates the location where the evaluation has been taken, the column indicates the 
#' function evaluated.
#' @description A Functional Object, represented respect to a Finite Element basis, is evaluated in a set of locations. 
#' @usage eval.FEM(FEM, locations)
#' @references Sangalli, L.M., Ramsay, J.O. & Ramsay, T.O., 2013. Spatial spline regression models. Journal of the Royal Statistical Society. Series B: Statistical Methodology, 75(4), pp.681.703. \cr
#'  Azzimonti, L. et al., 2014. Blood flow velocity field estimation via spatial regression with PDE penalization Blood flow velocity field estimation via spatial regression with PDE penalization. , (September), pp.37.41. \cr
#'  Devillers, O. et al. 2001. Walking in a Triangulation, Proceedings of the Seventeenth Annual Symposium on Computational Geometry

eval.FEM <- function(FEM, locations, CPP_CODE = TRUE)
{
  res = NULL
  if(CPP_CODE == FALSE)
  {
    res = R_eval.FEM(FEM, locations)
  }else
  {
    res = CPP_eval.FEM(FEM, locations, TRUE)
  }
  
  return(res)
}