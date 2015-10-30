#' Meuse river data set
#'
#' This data set gives locations and top soil heavy metal concentrations (ppm) of *****. 
#' For details on this dataset, see \code{meuse.all} in \url{https://cran.r-project.org/web/packages/gstat/gstat.pdf}. This data set is here provide, together with a definition of the domain boundary in {\link{MeuseBorder}}, as it is used in many examples to illustrate FEMr features. 
#'
#' \itemize{
#'   \item x. numeric vector indicating the x-coordinate.
#'   \item y. numeric vector indicating the y-coordinate.
#'   \item zinc. topsoil zinc concentration, ppm.
#'   ...
#' }
#'
#' @format A data frame with 155 rows and 12 variables.
#' @source \url{https://cran.r-project.org/web/packages/gstat/gstat.pdf}
#' @name MeuseData
NULL

#' Boundary of the Meuse River data set
#'
#' This file provides the boundary of the domain of the Meuse dataset.
#'
#' \itemize{
#'   \item V1. A vector having as entries the indices of the locations in {\link{MeuseData}} where a boundary segment starts from.
#'   \item V2. A vector having as entries the indices of the locations in {\link{MeuseData}} where a boundary segment ends to
#' }
#'
#' @format A data frame with 52 rows and 2 variables.
#' @name MeuseBorder
NULL

#' Simple mesh
#'
#' This contains a simple mesh, namely a TRIMESH2D object created with \code{create.MESH.2D}.
#'
#' @name mesh.2D.simple
NULL

#' Simple Rectangular mesh
#'
#' This contains a rectangular mesh, namely a TRIMESH2D object created with \code{create.MESH.2D}.
#'
#' @name mesh.2D.rectangular
NULL
