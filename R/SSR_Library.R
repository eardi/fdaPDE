#source("/home/eardi/workspace/SSR/R/SSR_RCode.R")
#source("/home/eardi/workspace/SSR/R/SSR_CPP.R")

create.FEM.basis = function(mesh, order, CPP_CODE = FALSE)
{
  #  CREATE.FEM.BASIS sets up a finite element basis for the analysis
  #  of spatial data.  It requires a triangular mesh as input.
  #  The triangular mesh is  the output of create.MESH.object or refine.MESH.object.
  #  The finite elements used for functional data analysis are first or second
  #  order Lagrangian elements.  These are triangles covering a region,
  #  and the basis system is piecewise polinomials (linear or quadratic). There is a basis
  #  function associated with each node in the system.
  #  When ORDER = 1 the basis system is piecewise linear and the nodes are the vertices of the triangles.
  #  When ORDER = 2 the basis system is piecewise quadratic and the nodes are points that are etiher 
  #  the vertices of the triangles or midpoints of edges of triangles.
  #
  #  Arguments:
  #  MESH..  A mesh object
  #  ORDER.  Order of elements, which may be either 1 or 2 (2 is default)
  #  
  #  Returns:
  #  An object of the basis class with parameters stored in member params,
  #  which in this case is a struct object with the mesh.
  type <- NULL
  if (class(mesh) == "TRIMESH2D")
  {
    type <- "FEM"
  }
  
  #  The number of basis functions corresponds to the number of vertices
  #  for order = 1, and to vertices plus edge midpoints for order = 2
  
  nbasis = dim(mesh$nodes)[[1]]
  
  eleProp = NULL
  if(CPP_CODE == FALSE)
  {
    eleProp = R_elementProperties(mesh, order)
  }
  
  basisobj = list(type = type, mesh = mesh, order = as.integer(order), nbasis = nbasis, J=eleProp$J, transf = eleProp$transf, metric = eleProp$metric)
  
  basisobj
}

smooth.FEM.basis<-function(locations = NULL, observations, basisobj, lambda, covariates = NULL, BC = NULL, GCV = TRUE, CPP_CODE = TRUE)
{
  bigsol = NULL  
  lambda = as.vector(lambda)
  
  if(CPP_CODE == FALSE)
  {
    print('R Code Execution')
    if(!is.null(BC))
    {
      print('Dirichlet Border Conditions will be ignored, set CPP_CODE = TRUE')
    }
    
    bigsol = R_smooth.FEM.basis(locations, observations, basisobj, lambda, covariates, GCV)   
  }else
  {
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.basis(locations, observations, basisobj, lambda, covariates, BC, GCV)
  }
  
  numnodes = nrow(basisobj$mesh$nodes)
  
  u = bigsol[[1]][1:numnodes,]
  s = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  felsplobj  = fobj(u, basisobj)
  laplacefd = fobj(s, basisobj)  
  
  reslist = NULL
  if(GCV == TRUE)
  {
    seq=getGCV(locations = locations, observations = observations, felsplobj = felsplobj, covariates = covariates, DOFseq = bigsol[[2]])
    reslist=list(felsplobj=felsplobj,laplacefd=laplacefd, DOF = bigsol[[2]], sigma2hatseq = seq$sigma2hatseq, GCVseq = seq$GCVseq)
  }else{
    reslist=list(felsplobj=felsplobj,laplacefd=laplacefd)
  }
  
  return(reslist)
}

smooth.FEM.PDE.basis<-function(locations = NULL, observations, basisobj, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = TRUE, CPP_CODE = TRUE)
{
  bigsol = NULL  
  lambda = as.vector(lambda)
  
  if(CPP_CODE == FALSE)
  {
    print('Function implemented only in C++, turn CPP_CODE = TRUE')  
  }else
  {
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.basis(locations, observations, basisobj, lambda, PDE_parameters, covariates, BC, GCV)
  }
  
  numnodes = nrow(basisobj$mesh$nodes)
  
  u = bigsol[[1]][1:numnodes,]
  s = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  felsplobj  = fobj(u, basisobj)
  laplacefd = fobj(s, basisobj)  
  
  reslist=list(felsplobj=felsplobj,laplacefd=laplacefd, DOF = bigsol[[2]])
  return(reslist)
}

smooth.FEM.PDE.SV.basis<-function(locations = NULL, observations, basisobj, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = TRUE, CPP_CODE = TRUE)
{
  bigsol = NULL  
  lambda = as.vector(lambda)
  
  if(CPP_CODE == FALSE)
  {
    print('Function implemented only in C++, turn CPP_CODE = TRUE')  
  }else
  {
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.SV.basis(locations, observations, basisobj, lambda, PDE_parameters, covariates, BC, GCV)
  }
  
  numnodes = nrow(basisobj$mesh$nodes)
  
  u = bigsol[[1]][1:numnodes,]
  s = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  felsplobj  = fobj(u, basisobj)
  laplacefd = fobj(s, basisobj)  
  
  reslist=list(felsplobj=felsplobj,laplacefd=laplacefd, DOF = bigsol[[2]])
  return(reslist)
}

fobj<-function(coefmat,basisobj)
{
  if(is.vector(coefmat))
  {
    coefmat = as.matrix(coefmat)
  }
  fclass = NULL
  fclass = list(coefmat=coefmat, basisobj=basisobj)
  class(fclass)<-"FOBJ"
  return(fclass)
}

eval.FEM.fobj <- function(fobj, locations, CPP_CODE = FALSE)
{
  # EVAL_FEM_FD evaluates the FEM fd object at points (X,Y)
  #
  #        arguments:
  # X         an array of x-coordinates.
  # Y         an array of y-coordinates.
  # FELSPLOBJ a FELspline object
  # FAST      a boolean indicating if the walking algorithm should be apply 
  #        output:
  # EVALMAT   an array of the same size as X and Y containing the value of 
  #           FELSPLOBJ at (X,Y).
  
  res = NULL
  if(CPP_CODE == FALSE)
  {
    res = R_eval.FEM.fobj(fobj, locations)
  }else
  {
    print("Not implemented yet")
  }
  
  return(res)
}

getGCV<-function(locations, observations, felsplobj, covariates = NULL, DOFseq)
{
  loc_nodes = NULL
  fnhat = NULL
  
  DOFseq = as.vector(DOFseq)
    
  if(is.null(locations))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    fnhat = felsplobj$coefmat[loc_nodes,]
  }else{
    loc_nodes = 1:length(observations)
    fnhat = eval.FEM.fobj(fobj = felsplobj, locations = locations, CPP_CODE = FALSE)
  }
  
  zhat = NULL
  zhat = matrix(nrow = length(loc_nodes), ncol = length(DOFseq))
  if(!is.null(covariates))
  {
    desmatprod = ( solve( t(covariates) %*% covariates ) ) %*% t(covariates)
    for ( i in 1:length(DOFseq))
    {
      betahat  = desmatprod %*% (observations-fnhat)[,i]
      zhat[,i]     = covariates %*% betahat + fnhat[,i]
    }
  }else{
    zhat = fnhat
  }
  
  np = length(loc_nodes)
  
  sigma2hatseq = numeric(length(DOFseq))
  GCVseq       = numeric(length(DOFseq))
  
  zhat <- as.matrix(zhat)
  
  for (i in 1:length(DOFseq))
  {
    sigma2hatseq[i] = t(observations[loc_nodes] - zhat[,i]) %*% (observations[loc_nodes] - zhat[,i]) / ( np - DOFseq[i] )
    GCVseq[i] = ( np / ( np - DOFseq[i] )) * sigma2hatseq[i]
  }
  
  return(list(sigma2hatseq = sigma2hatseq, GCVseq = GCVseq))
}
