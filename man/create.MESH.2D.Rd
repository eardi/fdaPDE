\name{create.MESH.2D}
\alias{create.MESH.2D}
\title{Creates a TRIMESH2D object}
\usage{
create.MESH.2D(nodes, nodesmarkers = NA, nodesattributes = NA, segments = NA, segmentsmarkers = NA, holes = NA, triangles = NA, order = 1, verbosity = 0)

}
\arguments{
  \item{nodes}{A n-by-2 matrix of nodes to be used as triangles' vertices.}

  \item{nodesmarkers}{A vector of length n specifying whether the node is on the border or not.}

  \item{nodesattributes}{A n-by-#attributes matrix specifying a set of attributes for each node.}
  
  \item{segments}{A matrix where each row is composed by two nodes indexes, specifying which segments need to be preserved in the triangulation. Commonly used for specifying the border.}
  
  \item{segmentsmarkers}{A matrix where the i-th row specifies a set of attributes of i-th segments.}
  
  \item{holes}{A matrix where the i-th row specifies the x and y coordinate of an hole.}
  
  \item{triangles}{A 3 or 6 column matrix, respectively for 1st and 2nd order meshes,
in which each row contains indices in \code{P} of
vertices. Useful to create a TRIMESH object when a triangulation has already been performed with an external software.}

  \item{order}{A numeric that can assume values \code{1,2}, that specifies the order of the elements of the traingulation}

  \item{verbosity}{A numeric that can assume values \code{0,1,2}, that specifies the output verbosity in the triangulation process.}
}
\value{
An object of class \code{TRIMESH2D}. This contains an updated list of the input parameters}

\description{
Creates a mesh object starting from a domain described by the input parameters. Whether only the nodes are specified the mesh is constructed on the convex hull. This function is based on the Triangle
library (\url{http://www.cs.cmu.edu/~quake/triangle.html}). The triangulation is a constrained conforming Delaunay triangulation in which additional vertices, called Steiner points, can be inserted into segments to improved the quality of the triangulation.
}

\examples{
## Creates an object TRIMESH2D on the convex hull of the specified nodes
mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)))
## Plot it
plot(mesh)

## Creates an object TRIMESH2D with a concavity and second order nodes
mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                     segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)))
## Plot it
plot(mesh)

## Creates an object TRIMESH2D with second order nodes starting from a first order triangulation
## specified by nodes and triangles
mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                     triangles=rbind(c(2,1,3), c(3,5,4), c(5,3,1)), order = 2)
## Plot it
plot(mesh)
}