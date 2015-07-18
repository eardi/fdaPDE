\name{smooth.FEM.PDE.basis}
\alias{smooth.FEM.PDE.basis}
\title{Compute a solution for a Spatial Regression with PDE Penalization}
\usage{
smooth.FEM.PDE.basis(locations = NULL, observations, basisobj, 
                    lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = TRUE, 
                    CPP_CODE = TRUE)
}
\arguments{

  \item{observations}{A vector specifying the observed values on the domain. The locations of the observations can be specified with the \code{locations} argument, otherwise the locations are intented to be the corresponding nodes of the mesh. \code{NA} values are admissible to indicate the missing value on the corresponding node.}
  
   \item{locations}{A 2 column matrix where each row specifies the coordinates of the corresponding observation.}

  \item{basisobj}{An onbject of type FEM; See
  \code{\link{create.FEM.basis}}.}
  
  \item{lambda}{A scalar smoothing parameter.}
  
  \item{PDE_parameters}{A list containing the parameters of the penalizing PDE, with: \code{K} a 2-by-2 matrix indicating the diffusion coefficient matrix, \code{beta} a vector of length 2 with the coefficients of the advection coefficients and \code{c} a numeric indicating the reaction coefficient.}
  
  \item{covariates}{A design matrix where each row represents the covariates associated to each row.}

  \item{BC}{A list with two vectors: \code{Indices}, a vector with the indices for the border points to apply a Dirichlet Border Condition;\code{Values} a vector with the values that the the nodes specified in \code{Indices} must assume.}

  \item{GCV}{If \code{TRUE} computes the trace of the smoothing matrix, the estimate of the error's variance and the Generalized Cross Validation parameter, for value of \code{lambda}.}
  
  \item{CPP_CODE}{if \code{TRUE} avoids the computation of some additional elements, not necessary if the functions working with the FEM basis are called with the flag \code{CPP_CODE=TRUE}}
}
\value{
A list with the following variables:
\item{\code{FELSPLOBJ}}{A FOBJ object of the FEM type defined by the coefficients vector resulting from smoothing.}
\item{\code{LAPLACEFD}}{A FOBJ object of the FEM type for the value of the Laplace operator}
\item{\code{DOF}}{If GCV is \code{TRUE}, a vector with the trace of the smoothing matrix for each lambda.}
\item{\code{sigma}}{If GCV is \code{TRUE}, a vector with the estimate of the standard deviation of the error for each lambda.}
\item{\code{GCV}}{If GCV is \code{TRUE}, a vector with the GCV index for each lambda.}

}

\description{
Compute a solution for a for a Spatial Regression with PDE Penalization model in Azzimonti et al (2013).
}

\examples{
library(RPDE)

data(MeuseData)
data(MeuseBorder)

order=1
mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = order)
plot(mesh)

data = log(MeuseData[,7])

basisobj = create.FEM.basis(mesh, order)

lambda = 10^3.5

PDE_parameters = list(K = 1*matrix(c(1,0,0,1), nrow = 2), beta = c(0,0), c = 0)

output = smooth.FEM.PDE.basis(observations = data, 
                              basisobj = basisobj, lambda = lambda, PDE_parameters = PDE_parameters)

plot(output$felsplobj)
}
