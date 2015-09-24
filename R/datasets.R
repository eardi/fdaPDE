#' Meuse river data set
#'
#' This data set gives locations and top soil heavy metal concentrations (ppm). It is used in many examples to illustrate FEMr features. 
#' For further details on this dataset, see \code{meuse.all} in \url{https://cran.r-project.org/web/packages/gstat/gstat.pdf}
#'
#' \itemize{
#'   \item x. numeric vector indicating the x-coordinate.
#'   \item y. numeric vector indicating the y-coordinate.
#'   \item zinc. topsoil zinc concentration, ppm.
#'   ...
#' }
#'
#' @format A data frame with 155 rows and 12 variables
#' @source \url{https://cran.r-project.org/web/packages/gstat/gstat.pdf}
#' @name MeuseData
NULL

#' Meuse Border description data set
#'
#' This file describes the boundaries of the domain of the Meuse dataset.
#'
#' \itemize{
#'   \item V1. An integer vector with the indices of the locations where the segments defining the border start from.
#'   \item V2. An integer vector with the indices of the locations where the segments defining the border end to
#' }
#'
#' @format A data frame with 52 rows and 2 variables
#' @name MeuseBorder
NULL

#' Simple mesh
#'
#' This contains a simple mesh, namely an object of the class TRIMESH2D, created with \code{create.MESH.2D}.
#'
#' @name mesh.2D.simple
NULL

#' Simple Rectangular mesh
#'
#' This contains a rectangular mesh, namely an object of the class TRIMESH2D, created with \code{create.MESH.2D}.
#'
#' @name mesh.2D.rectangular
NULL