triangulate_native <- function(P, PB, PA, S, SB,H, TR, flags) {
  ## It is necessary to check for NAs and NaNs, as the triangulate C
  ## code crashes if fed with them
  
  P  <- as.matrix(P)
  PB <- as.integer(PB)
  PA <- as.matrix(PA)
  S  <- as.matrix(S)
  SB <- as.integer(SB)
  H  <- as.matrix(H)
  TR  <- as.matrix(TR)
  
  storage.mode(P)  <- "double"
  storage.mode(PA) <- "double"
  storage.mode(PB) <- "integer"
  storage.mode(S)  <- "integer"
  storage.mode(SB) <- "integer"
  storage.mode(H)  <- "double"
  storage.mode(TR) <- "integer"
  storage.mode(flags) <- 'character'
  ## Call the main routine
  out <- .Call("R_triangulate_native",
               t(P),
               PB,
               PA,
               t(S),
               SB,
               H,
               t(TR),
               flags,
               PACKAGE="FEMr")
  names(out) <- c("P", "PB", "PA", "T", "S", "SB", "E", "EB","TN", "VP", "VE", "VN", "VA")
  class(out) <- "triangulation"
  return(out)
}

#' Create a Constrained Delaunay triangulation
#' 
#' @param nodes A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes
#' @param nodesmarkers A #nodes-vector containing '1' and '0' to indicate whether or not the node is a boundary node
#' @param nodesattributes A #nodes-by-n1 matrix containg some attributes to each point. 
#' These are copied unchanged to the output mesh for points in \code{nodes}. However for new points introducted by triangulation
#' procedure or further refinements, each new point added to the mesh will have quantities assigned to it by linear interpolation. This option is commonly
#' used to spread the Boundary Conditions to the new nodes introduced by the triangulation process
#' @param segments A #segments-by-2 matrix. Each row contains the indices (starting from zero) of the point where a segments start from and finisces to.
#' Segments are edges that persist after the triangulation. In he basic usage of the library segments are used to define the boundary
#' of the domain.
#' @param segmentsmarkers A #segments-vector containing '1' and '0' to indicate whether or not the segment is a boundary segment.
#' @param holes A #holes-by-2 matrix containg a point internal to each hole of the mesh. These points are used to carve the holes
#' after the triangulation procedure.
#' @param triangles A #triangles-by-3 or #triangles-by-6 matrix. This defines the triangles of the mesh already available. 
#' This option is usually used when a trianulation is already availble and is represented by the matrices \code{nodes} and \code{triangles}. 
#' However the \code{create.MESH.2D} should be used to produce a complete TRIMESH2D object. In https://www.cs.cmu.edu/~quake/triangle.highorder.html
#' a picture of the node ordering can be found.
#' @param order This can be '1' or '2'. It specifies if the triangular output elements should be represented by a 3 or 6 nodes. They are
#' respectively used for locally 1st and 2nd order Finite Elements.
#' @param verbosity This can be '0', '1' or '2'. It indicates the level of verbosity in triangulation process.
#' @description This function is a wrapper of the Triangle library (http://www.cs.cmu.edu/~quake/triangle.html). The function can be used
#' to create a triangulation starting from a list of points, to be used as triangles' vertices, and a list of segments defining the shape. However further options are availabe for more complex meshes. The resulting
#' triangulation is called Constrained Delaunay, thus is constructed in such a way to preserve the input \code{segments} without splitting them. 
#' @usage create.MESH.2D(nodes = nodeslist, segments = segmentlist)
#' @return An object of the class TRIMESH with the following variables:
#' \item{\code{nodes}}{A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes.}
#' \item{\code{nodesmarkers}}{A #nodes-vector containing '1' and '0' to indicate whether or not the node is a boundary node. This argument is rarely used.}
#' \item{\code{nodesattributes}}{A #nodes-by-n1 matrix containg some real attributes associated to each point.
#' These are copied unchanged to the output mesh for points in \code{nodes}. However for new points introducted by triangulation
#' procedure or further refinements, each new point added to the mesh will have quantities assigned to it by linear interpolation. This option is commonly
#' used to spread the Boundary Conditions to the new nodes introduced by the triangulation process.}
#' \item{\code{triangles}}{A #triangles-by-3 or #triangles-by-6 matrix, respectively where a 1st and 2nd order mesh is created. In https://www.cs.cmu.edu/~quake/triangle.highorder.html
#' a picture of the node ordering can be found.}
#' \item{\code{segments}}{A #segments-by-2 matrix. Each row contains the indices (starting from zero) of the point where a segments start from and finisces to.
#' Segments are edges that persist after the triangulation. In the basic usage of the library segments are used to define the boundary
#' of the domain.}
#' \item{\code{segmentsmarker}}{A #segments-vector containing '1' and '0' to indicate whether or not the segment is a boundary segment. This argument is rarely used.}
#' \item{\code{edges}}{A #edges-by-2 matrix. Each row contains the indices (starting from zero) of the point where an edge start from and finisces to.}
#' \item{\code{edgesmarkers}}{A #edges-vector containing '1' and '0' to indicate whether or not the segment is a boundary edge. This argument is rarely used.}
#' \item{\code{neighbors}}{A #triangles-by-3 matrix. Each row contains the indices of the three neighbouring triangles. '-1' if 
#' one side of the triangle is an edge on the boundary of the mesh.}
#' \item{\code{holes}}{A #holes-by-2 matrix containg a point internal to each hole of the mesh. These points are used to carve the holes
#' after the triangulation procedure.}
#' \item{\code{order}}{This can be '1' or '2'. It specifies if the triangular output elements is represented by a 3 or 6 nodes.}
#' @examples 
#' library(FEMr)
#' 
#' data(MeuseData)
#' data(MeuseBorder)
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = 1)
#' plot(mesh)
create.MESH.2D <- function(nodes, nodesmarkers = NA, nodesattributes = NA, segments = NA, segmentsmarkers = NA, holes = NA, triangles = NA, order = 1, verbosity = 0)
{ 
  ##########################
  ###   Input checking   ###
  ##########################
  
  ## Deal with nodes
  if (ncol(nodes) != 2) {
    stop("Matrix of nodes should have 2 columns")
  }
  
  ## Check that there are no duplicate rows in nodes
  if (anyDuplicated(nodes)) {
    stop("Duplicated nodes")
  }
  
  ## If attributes not specified, set them to a matrix with zero columns
  if (any(is.na(nodesattributes))) {
    nodesattributes <- matrix(0, nrow(nodes), 0)
  }
  ## Make sure the size of the point attribute matrix is correct
  if (nrow(nodesattributes) != nrow(nodes)) {
    stop("Point attribute matrix \'nodesattributes\' does not have same number of rows the point matrix \'nodes\'")
  }
  
  ## If boundary nodes not specified, set them to 0
  if (is.na(nodesmarkers)) {
    nodesmarkers <- 0
  }
  nodesmarkers <- rep(nodesmarkers, length.out=nrow(nodes))
  
  ## Deal with segments
  if (any(is.na(segments))) {
    segments <- matrix(0, 0, 2)
  } else {
    if (ncol(segments) != 2) {
      stop("Matrix of segments segments should have 2 columns")
    }
  }
  
  ## If boundary segments not specified, set them to 0
  if (any(is.na(segmentsmarkers))) {
    segmentsmarkers <- 0
  }
  segmentsmarkers <- rep(segmentsmarkers, length.out=nrow(segments))
  
  ## If hole not specified, set it to empty matrix
  if (any(is.na(holes))) {
    holes <- matrix(0, 0, 2)
  }
  
  ## If triangles are not already specified
  if(any(is.na(triangles))){
    triangles = matrix(0,nrow = 0, ncol = 3)
  }
  
  flags="ven" 
  
  if(nrow(segments) == 0){
    flags = paste(flags,"c",sep = '')
  }
  
  if(nrow(segments)>0){
    flags = paste(flags,"p",sep = '')
  }
  
  #If order=2 add flag for second order nodes
  if(order == 2){
    flags = paste(flags,"o2",sep = '')
  }
  if(order < 1 || order >2){
    print('Order must be 1 or 2')
  }
  
  if(nrow(triangles) > 0){
    flags = paste(flags,"r",sep = '')
  }
  
  if (verbosity == 0) {
    flags = paste(flags,"Q",sep = '')
  }
  if (verbosity == 1) {
    flags = paste(flags,"V",sep = '')
  }
  if (verbosity == 2) {
    flags = paste(flags,"VV",sep = '')
  }
  
  out<-NULL
  #If triangles is null it makes the trianglulation
  #If triangle is not null it makes a refinement with no parameter, to compose the mesh object
  out <- triangulate_native(  
    nodes,
    nodesmarkers,
    nodesattributes,
    segments,
    segmentsmarkers,
    holes,
    triangles,
    flags
  )
  
  names(out)[1]<-"nodes"
  names(out)[2]<-"nodesmarkers"
  names(out)[3]<-"nodesattributes"
  names(out)[4]<-"triangles"
  names(out)[5]<-"segments"
  names(out)[6]<-"segmentsmarkers"
  names(out)[7]<-"edges"
  names(out)[8]<-"edgesmarkers"
  names(out)[9]<-"neighbors"
  
  out[13]<-NULL
  out[12]<-NULL
  out[11]<-NULL
  out[10]<-NULL
  
  out[[10]] = holes
  names(out)[10]<-"holes"
  out[[11]] = order
  names(out)[11]<-"order"
  
  class(out)<-"TRIMESH2D"
  
  return(out)
}

#' Refine a Constrained Delaunay triangulation to a Conforming Delaunay triangulation
#' 
#' @param mesh An object of the class TRIMESH2D constructed through the function \code{create.MESH.2D}
#' @param minimum_angle A condition on the minimum angle that each vertex of each triangle of the output mesh should respect.
#' @param maximum_area A condition on the maximum area that each traingle should respect.
#' @param delaunay A boolean parameter indicating whether or not the output mesh should respect the Delaunay condition.
#' @param verbosity This can be '0', '1' or '2'. It indicates the level of verbosity in triangulation process.
#' @description This function is a wrapper of the Triangle library (http://www.cs.cmu.edu/~quake/triangle.html). It can be used to 
#' refine a mesh created previously with \link{create.MESH.2D}. The algorithm can add Steiner points (points through which the \code{segments} are splitted)
#' in order to meet the imposed conditions.
#' @usage refine.MESH.2D(mesh, minimum_angle, maximum_area, delaunay, verbosity)
#' @return An object of the class TRIMESH with the following variables:
#'  \item{\code{nodes}}{A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes.}
#'  \item{\code{nodesmarkers}}{A #nodes-vector containing '1' and '0' to indicate whether or not the node is a boundary node.}
#'  \item{\code{nodesattributes}}{A #nodes-by-n1 matrix containg some attributes for each point. 
#'        These are copied unchanged to the output mesh for points in \code{nodes}. However for new points introducted by triangulation
#'        procedure or further refinements, each new point added to the mesh will have quantities assigned to it by linear interpolation. This option is commonly
#'        used to spread the Boundary Conditions to the new nodes introduced by the triangulation process.}
#'  \item{\code{triangles}}{A #triangles-by-3 or #triangles-by-6 matrix, respectively where a 1st and 2nd order mesh is created. In https://www.cs.cmu.edu/~quake/triangle.highorder.html
#' a picture of the node ordering can be found.}
#'  \item{\code{segments}}{A #segments-by-2 matrix. Each row contains the indices (starting from zero) of the point where a segments start from and finisces to.
#'        Segments are edges that persist after the triangulation. In he basic usage of the library segments are used to define the boundary
#'        of the domain.}
#'  \item{\code{segmentsmarker}}{A #segments-vector containing '1' and '0' to indicate whether or not the segment is a boundary segment.}
#'  \item{\code{edges}}{A #edges-by-2 matrix. Each row contains the indices (starting from zero) of the point where an edge start from and finisces to.}
#'  \item{\code{edgesmarkers}}{A #edges-vector containing '1' and '0' to indicate whether or not the segment is a boundary edge.}
#'  \item{\code{neighbors}}{A #triangles-by-3 matrix. Each row contains the indices of the three neighbouring triangles. '-1' if 
#'  one side of the triangle is an edge on the boundary of the mesh.}
#'  \item{\code{holes}}{A #holes-by-2 matrix containg a point internal to each hole of the mesh. These points are used to carve the holes
#'        after the triangulation procedure.}
#'  
#'  \item{\code{order}}{This can be '1' or '2'. It specifies if the triangular output elements is represented by a 3 or 6 nodes.}
#' @examples 
#' library(FEMr)
#' 
#' data(MeuseData)
#' data(MeuseBorder)
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = 1)
#' plot(mesh)
#' 
#' mesh_refine <- refine.MESH.2D(mesh, minimum_angle = 30, maximum_area = 10000)
#' plot(mesh_refine)

refine.MESH.2D<-function(mesh, minimum_angle = NA, maximum_area = NA, delaunay = FALSE, verbosity = 0)
{ 
  flags="rpven" 
  
  if(!is.na(minimum_angle)){
    flags <- paste(flags, "q", minimum_angle, sep='')
  }
  
  if(!is.na(maximum_area)){
    flags <- paste(flags, "a", maximum_area, sep='')
  }
  
  if(delaunay){
    flags <- paste(flags, "D", sep='')
  }
  
  if(mesh$order==2){
    flags <- paste(flags, "o2", sep='')
  }
  
  if (verbosity == 0) {
    flags = paste(flags,"Q",sep = '')
  }
  if (verbosity == 1) {
    flags = paste(flags,"V",sep = '')
  }
  if (verbosity == 2) {
    flags = paste(flags,"VV",sep = '')
  }
  
  out<-NULL
  #If triangles is null it makes the trianglulation
  #If triangle is not null it makes a refinement with no parameter, to compose the mesh object
  out <- triangulate_native(
    mesh$nodes,
    mesh$nodesmarker,
    mesh$nodesattribute,
    mesh$segments,
    mesh$segmentmarker,
    mesh$holes,
    mesh$triangles,
    flags
  )
  
  names(out)[1]<-"nodes"
  names(out)[2]<-"nodesmarkers"
  names(out)[3]<-"nodesattributes"
  names(out)[4]<-"triangles"
  names(out)[5]<-"segments"
  names(out)[6]<-"segmentsmarkers"
  names(out)[7]<-"edges"
  names(out)[8]<-"edgesmarkers"
  names(out)[9]<-"neighbors"
  
  out[13]<-NULL
  out[12]<-NULL
  out[11]<-NULL
  out[10]<-NULL
  
  out[[10]] = mesh$holes
  names(out)[10]<-"holes"
  out[[11]] = order
  names(out)[11]<-"order"
  
  class(out)<-"TRIMESH2D"
  
  return(out)
}
#' Plot a mesh object of the class TRIMESH2D
#' 
#' @param x An object of the class TRIMESH2D defining the mesh. It can be created with \code{create.Mesh.2D} or \code{refine.Mesh.2D}.
#' @param ... Arguments representing the graphical parameters to be passed to methods, see \link[graphics]{par}.
#' @description Mesh objects created with the functions \code{create.MESH.2D} or \code{refine.MESH.2D} can be plotted.
#' @usage \method{plot}{TRIMESH2D}(x, ...)
#' @examples 
#' library(FEMr)
#' 
#' data(MeuseData)
#' data(MeuseBorder)
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = 1)
#' plot(mesh)
plot.TRIMESH2D<-function(x, ...)
{
  plot(x$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
  segments(x$nodes[x$edges[,1],1], x$nodes[x$edges[,1],2],
           x$nodes[x$edges[,2],1], x$nodes[x$edges[,2],2], ...)
  segments(x$nodes[x$segments[,1],1], x$nodes[x$segments[,1],2],
           x$nodes[x$segments[,2],1], x$nodes[x$segments[,2],2], col="red", ...)
}