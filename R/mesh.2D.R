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
  names(out)[2]<-"nodesmarker"
  names(out)[3]<-"nodesattribute"
  names(out)[4]<-"triangles"
  names(out)[5]<-"segments"
  names(out)[6]<-"segmentmarker"
  names(out)[7]<-"edges"
  names(out)[8]<-"edgemarker"
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
  names(out)[2]<-"nodesmarker"
  names(out)[3]<-"nodesattribute"
  names(out)[4]<-"triangles"
  names(out)[5]<-"segments"
  names(out)[6]<-"segmentmarker"
  names(out)[7]<-"edges"
  names(out)[8]<-"edgemarker"
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

plot.TRIMESH2D<-function(mesh, ...)
{
  plot(mesh$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
  segments(mesh$nodes[mesh$edges[,1],1], mesh$nodes[mesh$edges[,1],2],
           mesh$nodes[mesh$edges[,2],1], mesh$nodes[mesh$edges[,2],2], ...)
  segments(mesh$nodes[mesh$segments[,1],1], mesh$nodes[mesh$segments[,1],2],
           mesh$nodes[mesh$segments[,2],1], mesh$nodes[mesh$segments[,2],2], col="red", ...)
}
