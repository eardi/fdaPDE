#' Evaluate a FEM object at a set of point locations
#' 
#' @param locations A 2-colums matrix with the spatial locations where the FEM object should be evaluated.
#' @param A FEM object to be evaluated.
#' @param CPP_CODE Boolean. If \code{TRUE} the computation relies on the C++ implementation of a Visibility Walk Algorithm (Devillers et al. 2001). This usually ensures a fast computation.
#' @return 
#' A matrix of numeric evaluations of the \code{FEM} object. Each row indicates the location where the evaluation has been taken, the column indicates the 
#' function evaluated.
#' @description It evaluates a FEM object the specified set of locations.  
#' @usage eval.FEM(FEM, locations)
#' @references 
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
