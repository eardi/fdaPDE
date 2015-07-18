R_elementProperties=function(mesh, order)
{
  #MAKENODES produces:
  #  a matrix NODES containing coordinates for all of the nodes to be used, 
  #  a matrix NODEINDEX defining which nodes correspond to each element.  
  #  If NORDER is 2, the midpoint of each edge is computed and added to 
  #  POINTS to obtain matrix NODES. 
  #  The row index of that midpoint is then added to the rows of TRIANGLES 
  #  containing that edge to define NODEINDEX.
  #  If NORDER is 1, nodes corresponds to vertices and NODEINDEX is 
  #  identical to TRIANGLES.
  #
  #  nvert:  number of vertices
  #  nele:   number of triangles or elements
  #
  #Input: POINTS is an nvert by 2 matrix containing the x and y
  #   coordinates of each of the nvert points in the right-hand rule mesh.
  #   POINTS = P' where P is the points matrix for pde
  #   The call can use P directly (see below).
  #   TRIANGLES is T(1:3,:)' where T is the triangle index matrix from pde.
  #   Vertices must be numbered in the counterclockwise direction.
  #   NORDER is the order of elements, and may be either 1 or 2 (default)
  #
  #Output: 
  #     NODES:  a numnodes*2 matrix whose i'th row contains
  #       the coordinates of the i'th nodal variable.
  #       Nodes for the second order element consist of vertices and 
  #       midpoints of edges, that is, 6 per triangle.
  #       The first NVER rows of NODES is POINTS, and the remainder are
  #       the edge midpoints.
  #       Nodes for the first order element consist of only vertices.
  #
  #     NODEINDEX:  for NORDER == 2, an nele*6 matrix whose i'th row
  #       contains the row numbers (in NODES) of the
  #       nodal variables defining the i'th finite 
  #       element.  If the i'th row of FMESH is [V1 V2 V3]
  #       then the i'th row of nodeindex is
  #       [V1 V(12) V2 V(23) V3 V(31)], where Vi is the
  #       row number of the i'th point and V(ij) is the 
  #       row number of the midpoint of the edge defined
  #       by the i'th and j'th points.
  #       If NORDER == 1, NODEINDEX is TRIANGLES.
  #
  #  Last modified 4 February 2011 by Laura Sangalli.
  
  #  Eardi deleted the 'add second order nodes feature as the mesh generation can generate them' 
  #  adapted code to RTriangle nodes ordering (see Triangle on the web)
  #
  #  The first rows of nodes are the vertices
  
  
  nele = dim(mesh$triangles)[[1]]
  
  J   = matrix(0,nele,1)      #  vector of jacobian values
  metric = array(0,c(nele,2,2))  #  3-d array of metric matrices
  transf = array(0,c(nele,2,2))
  
  if (order ==1 || order == 2)
  {  
    for (i in 1:nele)
    {
      diff1x = mesh$nodes[mesh$triangles[i,2],1] - mesh$nodes[mesh$triangles[i,1],1]
      diff1y = mesh$nodes[mesh$triangles[i,2],2] - mesh$nodes[mesh$triangles[i,1],2]
      diff2x = mesh$nodes[mesh$triangles[i,3],1] - mesh$nodes[mesh$triangles[i,1],1]
      diff2y = mesh$nodes[mesh$triangles[i,3],2] - mesh$nodes[mesh$triangles[i,1],2]
      
      transf[i,,] = rbind(c(diff1x,diff2x),c(diff1y,diff2y))
      #  Jacobian or area of triangle
      J[i] = diff1x*diff2y - diff2x*diff1y
      
      #  Compute controvariant transformation matrix OSS: This is J^(-T)
      Ael = matrix(c(diff2y, -diff1y, -diff2x,  diff1x),nrow=2,ncol=2,byrow=T)/J[i]
      
      #  Compute metric matrix
      metric[i,,] = t(Ael)%*%Ael
    } 
  }else{
    stop("ORDER not 1 or 2")
  }
  
  FEStruct <- list(J=J, metric=metric, transf=transf)
  return(FEStruct)
}

R_mass=function(basisobj)
{
  #MASS produces the mass matrix containing integrals of products of
  #  nodal functions.  
  #
  #Input: NODESTRUCT is a struct object produced by function makenodes.
  #    It contains:
  #        ORDER     ... The order of the element (1 or 2)
  #        NODES     ... Coordinates of node points
  #        NODEINDEX ... indices of node points for each element
  #        J      ... Jacobian of the affine transformation of each
  #                      element to the master element
  #
  #Output: K0: the NNOD by NNOD matrix of sums of products of nodal basis
  #        functions.
  #        For each element i, the integral of the product 
  #        of the j'th and k'th shape functions over the i'th element is
  #        computed.  Then that value is the 
  #        (NODEINDEX(i,j),NODEINDEX(i,k))'th entry of the i'th elemental 
  #        mass matrix.
  
  
  nodes = basisobj$mesh$nodes
  triangles = basisobj$mesh$triangles
  J = basisobj$J
  order = basisobj$order
  
  nele  = dim(triangles)[1]
  nnod  = dim(nodes)[1]
  
  if (order < 1 || order > 2){
    stop("ORDER not 1 or 2")
  }
  
  K0M = NULL
  
  if (order ==2)
  {   
    #  the integrals of products of basis functions for master element:
    
    K0M = matrix(c( 6, -1, -1, -4,  0,  0,
                    -1,  6, -1,  0, -4,  0,
                    -1, -1,  6,  0,  0, -4,
                    -4,  0,  0, 32, 16, 16,
                    0, -4,  0, 16, 32, 16,
                    0,  0, -4, 16, 16, 32), ncol=6, nrow=6, byrow=T)/360
  }else if (order == 1)
  {  
    #  the integrals of products of basis functions for master element:
    K0M = matrix( c( 2,  1,  1,
                     1,  2,  1,
                     1 , 1,  2), ncol=3, nrow=3, byrow=T) / 24
    
    
  }
  
  # assemble the mass matrix
  K0 = matrix(0,nrow=nnod,ncol=nnod)
  for (el in 1:nele)
  {  
    ind = triangles[el,]
    K0[ind,ind] = K0[ind,ind] + K0M * J[el]
  }
  
  K0
}

R_stiff= function(basisobj)
{
  #STIFF1 produces the nnod*nnod stiffness matrix K1
  #defined (K1)jk = int(dpsik/da*dpsij/da + dpsik/db*dpsij/db).
  #
  #Input: NODESTRUCT is a struct object produced by function makenodes.
  #    It contains:
  #        ORDER     ... The order of the element (1 or 2)
  #        NODES     ... Coordinates of node points
  #        NODEINDEX ... indices of node points for each element
  #        J      ... Jacobian of the affine transformation of each
  #                      element to the master element
  #        METRIC    ... The crossproduct of the inverse of the linear
  #                      part of the transformation
  #
  #Output: K1 is an nnod*nnod matrix out which is
  #        the sum of the nele element stiffness matrices
  #        and the penalty stiffness matrix.
  #        These i'th element matrix has (ij)'th element defined
  #        as follows:  
  #        Let psita and psitb be the partial derivatives of the
  #        t'th shape function with respect to a and b (1<=t<=6).
  #        Then the integral of the sum of products
  #        (psija*psika+psijb+psikb) over the i'th element is
  #        computed.  Then that value is assigned to the
  #        (nodeindex(i,j),nodeindex(i,k))'th entry of the i'th elemental 
  #        stiffness matrix and the other elements are given the value zero.
  #
  #
  #  Last modified 4 February 2011 by Laura Sangalli.
  
  #  retrieve arrays from nodeStruct
  
  nodes = basisobj$mesh$nodes
  triangles = basisobj$mesh$triangles
  
  nele  = dim(triangles)[[1]]
  nnod  = dim(nodes)[[1]]
  J     = basisobj$J
  order = basisobj$order
  metric = basisobj$metric
  
  
  
  KXX = NULL
  KXY = NULL
  KYY = NULL
  
  if (order < 1 || order > 2){
    stop("ORDER not 1 or 2")
  }
  
  if (order == 2)
  {   
    #  values of K1 for master elements
    KXX = matrix( c( 3, 1, 0, 0, 0,-4, 
                     1, 3, 0, 0, 0,-4,
                     0, 0, 0, 0, 0, 0,
                     0, 0, 0, 8,-8, 0,
                     0, 0, 0,-8, 8, 0,
                     -4,-4, 0, 0, 0, 8), ncol=6, nrow=6, byrow=T)/6
    
    KXY = matrix(c( 3, 0, 1, 0,-4, 0, 
                    1, 0,-1, 4, 0,-4, 
                    0, 0, 0, 0, 0, 0,
                    0, 0, 4, 4,-4,-4,
                    0, 0,-4,-4, 4, 4,
                    -4, 0, 0,-4, 4, 4), ncol=6, nrow=6, byrow=T)/6
    
    KYY = matrix( c( 3, 0, 1, 0,-4, 0,
                     0, 0, 0, 0, 0, 0,
                     1, 0, 3, 0,-4, 0,
                     0, 0, 0, 8, 0,-8,
                     -4, 0,-4, 0, 8, 0,
                     0, 0, 0,-8, 0, 8), ncol=6, nrow=6, byrow=T)/6
  }else if (order == 1)
  {   
    KXX = matrix( c(  1, -1,  0,
                      -1,  1,  0,
                      0,  0,  0), ncol=3, nrow=3, byrow=T) /2
    
    KXY = matrix( c(  1,  0, -1,
                      -1,  0,  1,
                      0,  0,  0), ncol=3, nrow=3, byrow=T) /2
    
    KYY = matrix( c(  1,  0, -1,
                      0,  0,  0,
                      -1,  0,  1), ncol=3, nrow=3, byrow=T) /2
  }
  
  #  assemble the stiffness matrix
  K1   = matrix(0,nrow=nnod,ncol=nnod)
  for (el in 1:nele)
  {
    ind    = triangles[el,]
    K1M = (metric[el,1,1]*KXX    + metric[el,1,2]*KXY +
             metric[el,2,1]*t(KXY) + metric[el,2,2]*KYY)
    K1[ind,ind] = K1[ind,ind] + K1M*J[el]
  }
  
  K1
}

R_smooth.FEM.basis = function(locations, observations, basisobj, lambda, covariates = NULL, GCV)
{
  # SMOOTH.FEM.FD Compute a solution for a Spatial Spline problem 
  #
  #     Arguments:
  # FELSPLOBJ a FELspline object.
  # LAMBDA    a scalar smoothing parameter.
  # DATA      a n-by-2 set of noisy observations of the surface values.  
  #           DATA(:,1) indexes the points at which the 
  #           values in DATA(:,2) were observed.
  #
  #     Output:
  # FELSPLOBJ  ...  A FD object of the FEM type defined by the coefficient
  #                 vector resulting from smoothing
  # LAPLACEFD  ...  A FD object of the FEM type for the value of the 
  #                 Laplace operator if order == 2, or empty if order == 1
  #
  #
  #  Last modified on 8 February 2011 by Laura Sangalli
  
  
  
  #  check arguments
  
  
  #  check data argument  
  
  numnodes = nrow(basisobj$mesh$nodes)
  
  #  Construct penalty matrix and 'b' vector for Ax=b.
  
  #  ---------------------------------------------------------------
  # construct mass matrix K0 
  #  ---------------------------------------------------------------
  
  K0 = R_mass(basisobj)
  
  #  ---------------------------------------------------------------
  # construct stiffness matrix K1
  #  ---------------------------------------------------------------
  
  K1 = R_stiff(basisobj)
  
  
  #  ---------------------------------------------------------------
  # construct the penalty matrix P with ones on diagonal at data points
  #  ---------------------------------------------------------------
  
  #penalty = numeric(numnodes)
  #penalty[data[,1]] = 1
  #P = diag(penalty,nrow=numnodes)
  
  basismat = NULL
  Lmat = matrix(0,nrow = numnodes, ncol = numnodes)
  H = NULL
  Q = NULL
  
  if(!is.null(locations))
  {
    basismat = R_eval.FEM.basis(basisobj, locations)
  } 
  
  if(!is.null(covariates))
  {
    if(!is.null(locations))
    {
      H = covariates %*% ( solve( t(covariates) %*% covariates ) ) %*% t(covariates)
    }else{
      loc_nodes = (1:length(observations))[!is.na(observations)]
      H = covariates[loc_nodes,] %*% ( solve( t(covariates[loc_nodes,]) %*% covariates[loc_nodes,] ) ) %*% t(covariates[loc_nodes,])
    }
    Q=matrix(0,dim(H)[1], dim(H)[2])
    Q = diag(1,dim(H)[1])- H
  }
  
  if(!is.null(locations) && is.null(covariates))
  {
    Lmat = t(basismat) %*% basismat
  }
  if(!is.null(locations) && !is.null(covariates))
  {
    Lmat = t(basismat) %*% Q %*% basismat
  }
  if(is.null(locations) && is.null(covariates))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    Lmat[loc_nodes,loc_nodes]=diag(1,length(observations[loc_nodes]))
  }
  if(is.null(locations) && !is.null(covariates))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    Lmat[loc_nodes,loc_nodes] = Q
  }
  
  #print(Lmat[1:20,1:20])
  
  
  #  ---------------------------------------------------------------
  # construct vector b for system Ax=b
  #  ---------------------------------------------------------------
  
  b = matrix(numeric(numnodes*2),ncol=1)
  if(!is.null(locations) && is.null(covariates))
  {
    b[1:numnodes,] = t(basismat) %*% observations
  }
  if(!is.null(locations) && !is.null(covariates))
  {
    b[1:numnodes,] = t(basismat) %*% Q %*% observations
  }
  if(is.null(locations) && is.null(covariates))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    b[loc_nodes,] = observations[loc_nodes]
  }
  if(is.null(locations) && !is.null(covariates))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    b[loc_nodes,] = Q %*% observations[loc_nodes]
  }
  #print(b)
  
  #  ---------------------------------------------------------------
  # construct matrix A for system Ax=b.
  #  ---------------------------------------------------------------
  
  bigsol = list(solutions = matrix(nrow = 2*numnodes, ncol = length(lambda)), dof = vector(mode="numeric",length = length(lambda)))
  for(i in 1: length(lambda))
  {
    A  = rbind(cbind(Lmat, -lambda[i]*K1), cbind(K1, K0))
    #print(A)
    # solve system
    bigsol[[1]][,i] = solve(A,b)
    if(GCV == TRUE)
    {
      S = solve(Lmat + lambda[i] * K1 %*% solve(K0) %*% K1)
      if(is.null(locations))
      {
        loc_nodes = (1:length(observations))[!is.na(observations)]
        if(is.null(covariates))
        {
          bigsol[[2]][i]  = sum(diag(S[loc_nodes,loc_nodes]))
        }else{
          bigsol[[2]][i]  = ncol(covariates)+sum(diag(S[loc_nodes,loc_nodes]%*%Q))
        }
      }else{
        if(is.null(covariates))
        {
          bigsol[[2]][i]  = sum(diag(basismat%*%S[1:numnodes,1:numnodes]%*%t(basismat)))
        }else{
          bigsol[[2]][i]  = ncol(covariates)+sum(diag(basismat%*%S[1:numnodes,1:numnodes]%*%t(basismat)%*%Q))
        }
      }
    }else{
      bigsol[[2]][i]  = NA
    }
  }
  return(bigsol)
}

R_eval.FEM.basis <- function(basisobj, locations, nderivs = matrix(0,1,2))
{ 
  if(length(nderivs) != 2)
  {
    stop('NDERIVS not of length 2.')
  }
  
  if(sum(nderivs)>2)
  {
    stop('Maximum derivative order is greater than two.')
  }
  
  N = nrow(locations)
  nbasis = basisobj$nbasis
  
  # Augment Xvec and Yvec by ones for computing barycentric coordinates
  Pgpts = cbind(matrix(1,N,1),locations[,1],locations[,2])
  
  # Get nodes and index
  
  mesh = basisobj$mesh
  
  nodes = mesh$nodes
  triangles = mesh$triangles
  
  order = basisobj$order
  #nodeindex = params$nodeindex
  J = basisobj$J
  
  # 1st, 2nd, 3rd vertices of triangles
  
  v1 = nodes[triangles[,1],]
  v2 = nodes[triangles[,2],]
  v3 = nodes[triangles[,3],]
  
  if(order !=2 && order != 1)
  {
    stop('ORDER is neither 1 or 2.')
  }
  
  # Denominator of change of coordinates chsange matrix
  
  modJ = basisobj$J
  ones3 = matrix(1,3,1)
  modJMat = modJ %*% t(ones3)
  
  M1 = cbind(v2[,1]*v3[,2] - v3[,1]*v2[,2], v2[,2] - v3[,2], v3[,1] - v2[,1])/(modJMat)
  M2 = cbind(v3[,1]*v1[,2] - v1[,1]*v3[,2], v3[,2] - v1[,2], v1[,1] - v3[,1])/(modJMat)
  M3 = cbind(v1[,1]*v2[,2] - v2[,1]*v1[,2], v1[,2] - v2[,2], v2[,1] - v1[,1])/(modJMat)
  
  ind = matrix(0,N,1)
  for(i in 1:N)
  {
    ind[i] = R_insideIndex(mesh, locations[i,])
  }
  
  evalmat = matrix(0,N,nbasis)
  
  for(i in 1:N)
  {
    indi = ind[i]
    
    if(!is.nan(indi))
    {
      baryc1 = (M1[indi,]*Pgpts[i,]) %*% ones3
      baryc2 = (M2[indi,]*Pgpts[i,]) %*% ones3
      baryc3 = (M3[indi,]*Pgpts[i,]) %*% ones3
      
      if(order == 2)
      {
        if(sum(nderivs) == 0)
        {
          evalmat[i,triangles[indi,1]] = 2* baryc1^2 - baryc1
          evalmat[i,triangles[indi,2]] = 2* baryc2^2 - baryc2
          evalmat[i,triangles[indi,3]] = 2* baryc3^2 - baryc3
          evalmat[i,triangles[indi,6]] = 4* baryc1 * baryc2
          evalmat[i,triangles[indi,4]] = 4* baryc2 * baryc3
          evalmat[i,triangles[indi,5]] = 4* baryc3 * baryc1
        }
        else if(nderivs[1] == 1 && nderivs[2] == 0)
        {
          evalmat[i,triangles[indi,1]] = (4* baryc1 - 1) * M1[indi,2]
          evalmat[i,triangles[indi,2]] = (4* baryc2 - 1) * M2[indi,2]
          evalmat[i,triangles[indi,3]] = (4* baryc3 - 1) * M3[indi,2]
          evalmat[i,triangles[indi,6]] = (4* baryc2 ) * M1[indi,2] + 4*baryc1 * M2[indi,2]
          evalmat[i,triangles[indi,4]] = (4* baryc3 ) * M2[indi,2] + 4*baryc2 * M3[indi,2]
          evalmat[i,triangles[indi,5]] = (4* baryc1 ) * M3[indi,2] + 4*baryc3 * M1[indi,2]
        }
        else if(nderivs[1] == 0 && nderivs[2] == 1)
        {
          evalmat[i,triangles[indi,1]] = (4*baryc1 - 1)*M1[indi,3]
          evalmat[i,triangles[indi,2]] = (4*baryc2 - 1)*M2[indi,3]
          evalmat[i,triangles[indi,3]] = (4*baryc3 - 1)*M3[indi,3]
          evalmat[i,triangles[indi,6]] = 4*baryc2*M1[indi,3] + 4*baryc1*M2[indi,3]
          evalmat[i,triangles[indi,4]] = 4*baryc3*M2[indi,3] + 4*baryc2*M3[indi,3]
          evalmat[i,triangles[indi,5]] = 4*baryc1*M3[indi,3] + 4*baryc3*M1[indi,3]
        }
        else if(nderivs[1] == 1 && nderivs[2] == 1)
        {
          evalmat[i,triangles[indi,1]] = 4*M1[indi,2]%*%M1[indi,3];
          evalmat[i,triangles[indi,2]] = 4*M2[indi,2]%*%M2[indi,3];
          evalmat[i,triangles[indi,3]] = 4*M3[indi,2]%*%M3[indi,3];
          evalmat[i,triangles[indi,6]] = 4*M2[indi,2]%*%M1[indi,3] + 4*M2[indi,3]%*%M1[indi,2];
          evalmat[i,triangles[indi,4]] = 4*M3[indi,2]%*%M2[indi,3] + 4*M3[indi,3]%*%M2[indi,2];
          evalmat[i,triangles[indi,5]] = 4*M1[indi,2]%*%M3[indi,3] + 4*M1[indi,3]%*%M3[indi,2];
        }
        else if(nderivs[1] == 2 && nderivs[2] == 0)
        {
          evalmat[i,triangles[indi,1]] = 4*M1[indi,2]%*%M1[indi,2];
          evalmat[i,triangles[indi,2]] = 4*M2[indi,2]%*%M2[indi,2];
          evalmat[i,triangles[indi,3]] = 4*M3[indi,2]%*%M3[indi,2];
          evalmat[i,triangles[indi,6]] = 8*M2[indi,2]%*%M1[indi,2];
          evalmat[i,triangles[indi,4]] = 8*M3[indi,2]%*%M2[indi,2];
          evalmat[i,triangles[indi,5]] = 8*M1[indi,2]%*%M3[indi,2];
        }
        else if(nderivs[1] == 0 && nderivs[2] == 2)
        {
          evalmat[i,triangles[indi,1]] = 4*M1[indi,3]%*%M1[indi,3];
          evalmat[i,triangles[indi,2]] = 4*M2[indi,3]%*%M2[indi,3];
          evalmat[i,triangles[indi,3]] = 4*M3[indi,3]%*%M3[indi,3];
          evalmat[i,triangles[indi,6]] = 8*M2[indi,3]%*%M1[indi,3];
          evalmat[i,triangles[indi,4]] = 8*M3[indi,3]%*%M2[indi,3];
          evalmat[i,triangles[indi,5]] = 8*M1[indi,3]%*%M3[indi,3];
        }
        else
        {
          stop('Inadmissible derivative orders.')
        }
      }
      else
      {
        if(sum(nderivs) == 0)
        {
          evalmat[i,triangles[indi,1]] = baryc1;
          evalmat[i,triangles[indi,2]] = baryc2;
          evalmat[i,triangles[indi,3]] = baryc3;
        }
        else if(nderivs[1] == 1 && nderivs[2] == 0)
        {
          evalmat[i,triangles[indi,1]] = M1[indi[1],2];
          evalmat[i,triangles[indi,2]] = M1[indi[2],2];
          evalmat[i,triangles[indi,3]] = M1[indi[3],2]; 
        }
        else if(nderivs[1] == 0 && nderivs[2] == 1)
        {
          evalmat[i,triangles[indi,1]] = M1[indi[1],3];
          evalmat[i,triangles[indi,2]] = M1[indi[2],3];
          evalmat[i,triangles[indi,3]] = M1[indi[3],3]; 
        }
        else
        {
          stop('Inadmissible derivative orders.')
        }
      }
    }
  }
  return(evalmat)
}

R_eval.FEM.fobj <- function(fobj, locations)
{ 
  if (is.vector(locations))
  {
    locations = t(as.matrix(locations))
  }
  
  N = nrow(locations)
  
  # Augment Xvec and Yvec by ones for computing barycentric coordinates
  Pgpts = cbind(matrix(1,N,1),locations[,1],locations[,2])
  
  # Get nodes and index
  
  mesh = basisobj$mesh
  
  nodes = mesh$nodes
  triangles = mesh$triangles
  coefmat = fobj$coefmat
  nsurf = dim(coefmat)[2]
  
  basisobj = fobj$basisobj
  order = basisobj$order
  #nodeindex = params$nodeindex
  J = basisobj$J
  
  # 1st, 2nd, 3rd vertices of triangles
  
  v1 = nodes[triangles[,1],]
  v2 = nodes[triangles[,2],]
  v3 = nodes[triangles[,3],]
  
  if(order !=2 && order != 1)
  {
    stop('ORDER is neither 1 or 2.')
  }
  
  # Denominator of change of coordinates chsange matrix
  
  modJ = basisobj$J
  ones3 = matrix(1,3,1)
  modJMat = modJ %*% t(ones3)
  
  M1 = cbind(v2[,1]*v3[,2] - v3[,1]*v2[,2], v2[,2] - v3[,2], v3[,1] - v2[,1])/(modJMat)
  M2 = cbind(v3[,1]*v1[,2] - v1[,1]*v3[,2], v3[,2] - v1[,2], v1[,1] - v3[,1])/(modJMat)
  M3 = cbind(v1[,1]*v2[,2] - v2[,1]*v1[,2], v1[,2] - v2[,2], v2[,1] - v1[,1])/(modJMat)
  
  ind = matrix(0,N,1)
  for(i in 1:N)
  {
    ind[i] = R_insideIndex(mesh, as.numeric(locations[i,]))
  }
  
  evalmat = matrix(NA, nrow=N, ncol=nsurf)
  
  for (isurf in 1:nsurf)
  {
    for(i in 1:N)
    {
      indi = ind[i]
      
      if(!is.nan(indi))
      {
        baryc1 = (M1[indi,]*Pgpts[i,]) %*% ones3
        baryc2 = (M2[indi,]*Pgpts[i,]) %*% ones3
        baryc3 = (M3[indi,]*Pgpts[i,]) %*% ones3
        
        if(order == 2)
        {
          c1 = coefmat[triangles[indi,1],isurf]
          c2 = coefmat[triangles[indi,2],isurf]
          c3 = coefmat[triangles[indi,3],isurf]
          c4 = coefmat[triangles[indi,6],isurf]
          c5 = coefmat[triangles[indi,4],isurf]
          c6 = coefmat[triangles[indi,5],isurf]
          
          fval =  c1*(2* baryc1^2 - baryc1) +
            c2*(2* baryc2^2 - baryc2) +
            c3*(2* baryc3^2 - baryc3) +
            c4*(4* baryc1 * baryc2) +
            c5*(4* baryc2 * baryc3) +
            c6*(4* baryc3 * baryc1)
          evalmat[i,isurf] = fval
        }
        else
        {
          c1 = coefmat[triangles[indi,1],isurf]
          c2 = coefmat[triangles[indi,2],isurf]
          c3 = coefmat[triangles[indi,3],isurf]
          fval = c1*baryc1 + c2*baryc2 + c3*baryc3
          evalmat[i,isurf] = fval
        }
      }
    }
  }
  return(evalmat)
}

R_tricoefCal <- function(mesh)
{
  UseMethod("R_tricoefCal",mesh)
}

R_tricoefCal.TRIMESH2D = function(mesh)
{
  #  TRICOEFCAL compute the coefficient matrix TRICOEF
  #  required to test of a point is indside a triangle
  
  #  Last modified 24 June 2010 by Laura Sangalli
  
  nodes = mesh$nodes
  triangles = mesh$triangles
  
  ntri   = dim(triangles)[[1]]
  
  #  compute coefficients for computing barycentric coordinates if
  #  needed
  
  tricoef = matrix(0, nrow=ntri, ncol=4)
  tricoef[,1] = nodes[triangles[,1],1]-nodes[triangles[,3],1]
  tricoef[,2] = nodes[triangles[,2],1]-nodes[triangles[,3],1]
  tricoef[,3] = nodes[triangles[,1],2]-nodes[triangles[,3],2]
  tricoef[,4] = nodes[triangles[,2],2]-nodes[triangles[,3],2]
  detT = matrix((tricoef[,1]*tricoef[,4] - tricoef[,2]*tricoef[,3]),ncol=1)
  tricoef = tricoef/(detT %*% matrix(1,nrow=1,ncol=4))
  
  return(tricoef)
}

R_tricoefCal.default = function(mesh)
{
  print("Mesh not supported in R")
}

R_insideIndex <- function(mesh, location)
{
  UseMethod("R_insideIndex",mesh)
}

R_insideIndex.TRIMESH2D = function (mesh, location)
{
  #  insideIndex returns the index of the triangle containing the point
  # (X,Y) if such a triangle exists, and NaN otherwise.
  #  TRICOEF may have already been calculated for efficiency,
  #  but if the function is called with four arguments, it is calculated.
  
  
  #  Last modified 24 June 2010 by Laura Sangalli
  
  eps=2.2204e-016
  small = 10000*eps
  
  nodes = mesh$nodes
  triangles = mesh$triangles
  X = location[1]
  Y = location[2]
  
  ntri   = dim(triangles)[[1]]
  indtri   = matrix(1:ntri,ncol=1)
  
  #  compute coefficients for computing barycentric coordinates if needed
  
  tricoef = R_tricoefCal(mesh)
  
  #  compute barycentric coordinates
  r3 = X - nodes[triangles[,3],1]
  s3 = Y - nodes[triangles[,3],2]
  lam1 = ( tricoef[,4]*r3 - tricoef[,2]*s3)
  lam2 = (-tricoef[,3]*r3 + tricoef[,1]*s3)
  lam3 = 1 - lam1 - lam2
  
  #  test these coordinates for a triple that are all between 0 and 1
  int  = (-small <= lam1 & lam1 <= 1+small) & 
    (-small <= lam2 & lam2 <= 1+small) & 
    (-small <= lam3 & lam3 <= 1+small)
  
  #  return the index of this triple, or NaN if it doesn't exist
  indi = indtri[int]
  if (length(indi)<1)
  {
    ind = NA
  }else{
    ind = min(indi)
  }
  
  ind
}

R_insideIndex.default = function(mesh, location)
{
  print("Mesh not supported in R")
}

R_eval_local.FEM.fobj = function(fobj, locations, element_index)
{
  N = nrow(locations)
  # Augment Xvec and Yvec by ones for computing barycentric coordinates
  Pgpts = cbind(matrix(1,N,1),locations[,1],locations[,2])
  
  # Get nodes and index
  
  mesh = basisobj$mesh
  nodes = mesh$nodes
  triangles = mesh$triangles
  coefmat = fobj$coefmat
  nsurf = dim(coefmat)[2]
  
  basisobj = fobj$basisobj
  order = basisobj$order
  #nodeindex = params$nodeindex
  J = basisobj$J
  
  # 1st, 2nd, 3rd vertices of triangles
  
  v1 = nodes[triangles[element_index,1],]
  v2 = nodes[triangles[element_index,2],]
  v3 = nodes[triangles[element_index,3],]
  
  if(order !=2 && order != 1)
  {
    stop('ORDER is neither 1 or 2.')
  }
  
  # Denominator of change of coordinates chsange matrix
  
  modJ = basisobj$J[element_index]
  ones3 = matrix(1,3,1)
  #modJMat = modJ %*% t(ones3)
  
  M1 = c(v2[1]*v3[2] - v3[1]*v2[2], v2[2] - v3[2], v3[1] - v2[1])/(modJ)
  M2 = c(v3[1]*v1[2] - v1[1]*v3[2], v3[2] - v1[2], v1[1] - v3[1])/(modJ)
  M3 = c(v1[1]*v2[2] - v2[1]*v1[2], v1[2] - v2[2], v2[1] - v1[1])/(modJ)
  
  evalmat = matrix(NA, nrow=N, ncol=nsurf)
  
  for (isurf in 1:nsurf)
  {
    for(i in 1:N)
    {
      baryc1 = (M1*Pgpts[i,]) %*% ones3
      baryc2 = (M2*Pgpts[i,]) %*% ones3
      baryc3 = (M3*Pgpts[i,]) %*% ones3
      
      if(order == 2)
      {
        c1 = coefmat[triangles[element_index,1],isurf]
        c2 = coefmat[triangles[element_index,2],isurf]
        c3 = coefmat[triangles[element_index,3],isurf]
        c4 = coefmat[triangles[element_index,6],isurf]
        c5 = coefmat[triangles[element_index,4],isurf]
        c6 = coefmat[triangles[element_index,5],isurf]
        
        fval =  c1*(2* baryc1^2 - baryc1) +
          c2*(2* baryc2^2 - baryc2) +
          c3*(2* baryc3^2 - baryc3) +
          c4*(4* baryc1 * baryc2) +
          c5*(4* baryc2 * baryc3) +
          c6*(4* baryc3 * baryc1)
        evalmat[i,isurf] = fval
      }else{
        c1 = coefmat[triangles[element_index,1],isurf]
        c2 = coefmat[triangles[element_index,2],isurf]
        c3 = coefmat[triangles[element_index,3],isurf]
        fval = c1*baryc1 + c2*baryc2 + c3*baryc3
        evalmat[i,isurf] = fval
      }
    }
  }
  return(evalmat)
}

R_plot.ORD1.FOBJ = function(fobj)  
{
  # PLOT  Plots a FEM object FDOBJ over a rectangular grid defined by 
  # vectors X and Y;
  #
  
  #  Last modified 4 February 2011 by Laura Sangalli.
  
  
  #   if (!is.fd(fdobj))
  #   {
  #     stop('FDOBJ is not an FD object')
  #   }
  
  nodes = fobj$basisobj$mesh$nodes
  triangles = fobj$basisobj$mesh$triangles
  
  coefmat = fobj$coefmat
  
  basisobj = fobj$basis
  
  mesh = basisobj$mesh
  
  heat = heat.colors(100)
  
  nsurf = dim(coefmat)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d()
    axes3d()
    
    z = coefmat[as.vector(t(triangles)),isurf]
    rgl.triangles(x = nodes[as.vector(t(triangles)) ,1], y = nodes[as.vector(t(triangles)) ,2], 
                  z=coefmat[as.vector(t(triangles)),isurf], 
                  color = heat[round(99*(z- min(z))/(max(z)-min(z)))+1])
    aspect3d(2,2,1)
    rgl.viewpoint(0,-45)
    if (nsurf > 1)
    {readline("Press a button for the next plot...")}
  }
}

R_plot.ORDN.FOBJ = function(fobj, num_refinements)  
{
  coefmat = fobj$coefmat
  
  basisobj = fobj$basis
  
  mesh = basisobj$mesh
  
  heat = heat.colors(100)
  
  coefmat = fobj$coefmat
  
  # num_refinements sets the number od division on each triangle edge to be applied for rifenment
  if(is.null(num_refinements))
  {
    num_refinements = 20
  }
  
  # For the reference triangles we construct a regular mesh
  x = seq(from = 0, to = 1, length.out = num_refinements+1)
  y = seq(from = 0, to = 1, length.out = num_refinements+1)
  points_ref = expand.grid(x,y)
  points_ref = points_ref[points_ref[,1] + points_ref[,2] <= 1,]
  
  
  meshi = create.MESH.2D(nodes = points_ref, order = 1)   
  #plot(meshi)
  
  # locations is the matrix with that will contain the coordinate of the points where the function is 
  # evaluated (1st and 2nd column) and the columns with the evaluation of the ith fucntion on that point
  
  locations = matrix(nrow = nrow(mesh$triangles)*nrow(meshi$nodes), ncol = 2+ncol(coefmat))
  triangles = matrix(nrow = nrow(mesh$triangles)*nrow(meshi$triangles), ncol = 3)
  tot = 0
  
  
  for (i in 1:nrow(mesh$triangles))
  {
    # For each traingle we define a fine mesh as the transofrmation of the one constructed for the reference
    pointsi = t(basisobj$transf[i,,]%*%t(meshi$nodes) + mesh$nodes[mesh$triangles[i,1],])
    #We evaluate the fine mesh OBS: we know the triangle we are working on no need for point location
    z = R_eval_local.FEM.fobj(fobj, locations = pointsi, element_index = i)
    
    #We store the results
    locations[((i-1)*nrow(pointsi)+1):(i*nrow(pointsi)),] = cbind(pointsi,z)
    triangles[((i-1)*nrow(meshi$triangles)+1):(i*nrow(meshi$triangles)),] = meshi$triangles+tot
    tot = tot + nrow(meshi$nodes)
  }
  
  heat = heat.colors(100)
  
  nsurf = dim(coefmat)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d()
    axes3d()
    rgl.pop("lights") 
    light3d(specular="black") 
    z = locations[as.vector(t(triangles)), 2 + isurf]
    rgl.triangles(x = locations[as.vector(t(triangles)) ,1], y = locations[as.vector(t(triangles)) ,2], 
                  z = z, 
                  color = heat[round(99*(z-min(z))/(max(z)-min(z)))+1])
    aspect3d(2,2,1)
    rgl.viewpoint(0,-45)
    if (nsurf > 1)
    {readline("Press a button for the next plot...")}
  }
}

R_plot_old.ORDN.FOBJ = function(fobj, num_refinements)  
{
  coefmat = fobj$coefmat
  
  basisobj = fobj$basis
  
  mesh = basisobj$mesh
  
  heat = heat.colors(100)
  
  coefmat = fobj$coefmat
  
  locations = NULL
  triangles = NULL
  tot = 0
  
  maximum_area = NULL
  if(is.null(num_refinements))
  {
    maximum_area = fobj$basisobj$J / 500
  }else{
    maximum_area = fobj$basisobj$J / (2*num_refinements)
  }
  
  for (i in 1:nrow(mesh$triangles))
  {
    trianglei= t(as.matrix(mesh$triangles[i,1:3]))
    pointsi = mesh$nodes[mesh$triangles[i,1:3],]
    
    meshi = create.MESH.2D(nodes = pointsi, order = 1)    
    meshi2 = refine.MESH.2D(meshi, maximum_area = maximum_area[i])
    z = R_eval_local.FEM.fobj(fobj, locations = meshi2$nodes, element_index = i)
    locations = rbind(locations, cbind(meshi2$nodes,z))
    
    triangles = rbind(triangles, meshi2$triangles+tot) 
    tot = tot + nrow(meshi2$nodes)
  }
  
  heat = heat.colors(100)
  
  nsurf = dim(coefmat)[[2]]
  for (isurf in 1:nsurf)
  {
    axes3d()
    rgl.triangles(x = locations[as.vector(t(triangles)) ,1], y = locations[as.vector(t(triangles)) ,2], 
                  z= locations[as.vector(t(triangles)), 2 + isurf], 
                  color = heat[round(99*(locations[as.vector(t(triangles)), 3]- min(locations[as.vector(t(triangles)), 3]))/(max(locations[as.vector(t(triangles)), 3])-min(locations[as.vector(t(triangles)), 3])))+1])
    aspect3d(2,2,1)
    rgl.viewpoint(0,-45)
    if (nsurf > 1)
    {readline("Press a button for the next plot...")}
  }
}

plot.FOBJ = function(fobj, num_refinements = NULL)  
{
  if(fobj$basisobj$order == 1)
  {
    R_plot.ORD1.FOBJ(fobj)
  }else{
    R_plot.ORDN.FOBJ(fobj, num_refinements)
  }
}

R_image.ORD1.FOBJ = function(fobj)  
{
  # PLOT  Plots a FEM object FDOBJ over a rectangular grid defined by 
  # vectors X and Y;
  #
  
  #  Last modified 4 February 2011 by Laura Sangalli.
  
  
  #   if (!is.fd(fdobj))
  #   {
  #     stop('FDOBJ is not an FD object')
  #   }
  
  nodes = fobj$basisobj$mesh$nodes
  triangles = fobj$basisobj$mesh$triangles
  
  coefmat = fobj$coefmat
  
  basisobj = fobj$basis
  
  mesh = basisobj$mesh
  
  heat = heat.colors(100)
  
  nsurf = dim(coefmat)[[2]]
  for (isurf in 1:nsurf)
  {
    #rgl.open()
    axes3d()
    rgl.pop("lights") 
    light3d(specular="black") 
    z = coefmat[as.vector(t(triangles)),isurf]
    rgl.triangles(x = nodes[as.vector(t(triangles)) ,1], y = nodes[as.vector(t(triangles)) ,2], 
                  z=0, 
                  color = heat[round(99*(z- min(z))/(max(z)-min(z)))+1])
    aspect3d(2,2,1)
    rgl.viewpoint(0,0)
    if (nsurf > 1)
    {readline("Press a button for the next plot...")}
  }
}

R_image.ORDN.FOBJ = function(fobj, num_refinements)  
{
  coefmat = fobj$coefmat
  
  basisobj = fobj$basis
  
  mesh = basisobj$mesh
  
  heat = heat.colors(100)
  
  coefmat = fobj$coefmat
  
  if(is.null(num_refinements))
  {
    num_refinements = 20
  }
  
  x = seq(from = 0, to = 1, length.out = num_refinements+1)
  y = seq(from = 0, to = 1, length.out = num_refinements+1)
  points_ref = expand.grid(x,y)
  points_ref = points_ref[points_ref[,1] + points_ref[,2] <= 1,]
  
  meshi = create.MESH.2D(nodes = points_ref, order = 1)   
  #plot(meshi)
  
  locations = matrix(nrow = nrow(mesh$triangles)*nrow(meshi$nodes), ncol = 3*ncol(coefmat))
  triangles = matrix(nrow = nrow(mesh$triangles)*nrow(meshi$triangles), ncol = 3*ncol(coefmat))
  tot = 0
  
  for (i in 1:nrow(mesh$triangles))
  {
    pointsi = t(basisobj$transf[i,,]%*%t(meshi$nodes) + mesh$nodes[mesh$triangles[i,1],])
    z = R_eval_local.FEM.fobj(fobj, locations = pointsi, element_index = i)
    
    #mesh3 <- addNormals(subdivision3d(tmesh3d(vertices = t(locations), indices = t(triangles), homogeneous = FALSE)),deform = TRUE)
    
    locations[((i-1)*nrow(pointsi)+1):(i*nrow(pointsi)),] = cbind(pointsi,z)
    triangles[((i-1)*nrow(meshi$triangles)+1):(i*nrow(meshi$triangles)),] = meshi$triangles+tot
    tot = tot + nrow(meshi$nodes)
  }
  
  heat = heat.colors(100)
  
  nsurf = dim(coefmat)[[2]]
  for (isurf in 1:nsurf)
  {
    axes3d()
    rgl.pop("lights") 
    light3d(specular="black") 
    z = locations[as.vector(t(triangles)), 2 + isurf];
    rgl.triangles(x = locations[as.vector(t(triangles)) ,1], y = locations[as.vector(t(triangles)) ,2], 
                  z=0, 
                  color = heat[round(99*(z- min(z))/(max(z)-min(z)))+1])
    aspect3d(2,2,1)
    rgl.viewpoint(0,0)
    if (nsurf > 1)
    {readline("Press a button for the next plot...")}
  }
}

image.FOBJ = function(fobj, num_refinements = NULL)  
{
  if(fobj$basisobj$order == 1)
  {
    R_image.ORD1.FOBJ(fobj)
  }else{
    R_image.ORDN.FOBJ(fobj, num_refinements)
  }
}