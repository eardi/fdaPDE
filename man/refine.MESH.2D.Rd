\name{refine.MESH.2D}
\alias{refine.MESH.2D}
\title{Refine the triangulation}
\usage{
refine.MESH.2D(mesh, minimum_angle = NA, maximum_area = NA, delaunay = FALSE, verbosity = 0)
}
\arguments{
  \item{mesh}{A TRIMESH2D object created through \code{\link{create.MESH.2D}}.}

  \item{minimum_angle}{A numeric specifying the minimum angle of the triangles in the ouput triangulation.}

  \item{maximum_area}{A numeric specifying the maximum area of the triangles in the ouput triangulation.}
  
  \item{delaunay}{If \code{TRUE} the output triangulation is a Delaunay triangulation}
  
  \item{verbosity}{A numeric that can assume values \code{0,1,2}, that specifies the output verbosity in the triangulation process.}
}
\value{
An object of class \code{TRIMESH2D}.}

\description{
Refines a TRIMESH2D object following the constrained imposed with the input parameters. This function is based on the Triangle library (\url{http://www.cs.cmu.edu/~quake/triangle.html}). The triangulation is a constrained conforming Delaunay triangulation in which additional vertices, called Steiner points, can be inserted into segments to improved the quality of the triangulation.
}

\examples{
## Creates an object TRIMESH2D with a concavity and second order nodes
mesh_coarse<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                     segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)))
## Plot it
plot(mesh_coarse)                   
## Refines the the triangulation in specified in the \code{mesh_coarse} object
mesh<-refine.MESH.2D(mesh_coarse,maximum_area = 0.005, delaunay = TRUE)
## Plot the refined mesh
plot(mesh)
}