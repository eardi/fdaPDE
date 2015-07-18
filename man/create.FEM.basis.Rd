\name{create.FEM.basis}
\alias{create.FEM.basis}
\title{Creates a Finite Element Method basis}
\usage{
create.FEM.basis(mesh, order, CPP_CODE = FALSE)
}
\arguments{

  \item{mesh}{A \code{TRIMESH2D} object; see
  \code{\link{create.MESH.2D}}.}

  \item{order}{Order of elements, which may be either 1 or 2. The order of the element must be less or equal respect to the order of the mesh.}

  \item{CPP_CODE}{if \code{TRUE} avoids the computation of some additional elements, not necessary if the functions working with the FEM basis are called with the flag \code{CPP_CODE=TRUE}}
}
\value{
A object of class \code{FEM}. This contains the
mesh along with the following variables:

\item{\code{order}}{Order of elements.} \item{\code{nbasis}}{The number of basis.} \item{\code{J}}{The area of each triangle of the basis.} \item{\code{transf}}{A matrix such that J[i,,] is the 2-by-2 tranformation matrix such that transforms the nodes of the reference triangle in the node of the i-th triangle.}
\item{\code{metric}}{A matrix such that metric[i,,] is the 2-by-2 tranformation matrix equal to \code{J[i,,]^{-1}*J[i,,]^{-T}}
}
}

\description{
Sets up a finite element basis for the analysis of spatial data. It requires a TRIMESH2D object as input.
The finite elements used for functional data analysis are first or second order Lagrangian elements. These are triangles covering a region and the basis system is piecewise polinomials (linear or quadratic). There is a basis  function associated with each node in the system.
When ORDER = 1 the basis system is piecewise linear and the nodes are the vertices of the triangles.
When ORDER = 2 the basis system is piecewise quadratic and the nodes are points that are etiher the vertices of the triangles or midpoints of edges of triangles.
}

\examples{
## Creates an object TRIMESH2D with a concavity and second order nodes
mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                     segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)), order=1)
## Plot it
plot(mesh)                   

# Creates the basis object
basisobj = create.FEM.basis(mesh, order=1)
}
