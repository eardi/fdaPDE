#' Compute a solution for a Spatial Spline problem
#' 
#' @param observations A vector specifying the observed values on the domain. 
#' The locations of the observations can be specified with the \code{locations} 
#' argument, otherwise the locations are intented to be the corresponding nodes of the mesh. 
#' \code{NA} values are admissible to indicate the missing value on the corresponding node.
#' @param locations A 2 column matrix where each row specifies the coordinates of the corresponding observation.
#' @param basisobj An an object of type FEM; See \code{\link{create.FEM.basis}}.
#' @param lambda A scalar smoothing parameter.
#' @param covariates A design matrix where each row represents the covariates associated to each row.
#' @param BC A list with two vectors: 
#'        \code{Indices}, a vector with the indices for the border points to apply a Dirichlet Border Condition;
#'        \code{Values} a vector with the values that the the nodes specified in \code{Indices} must assume.
#' @param GCV If \code{TRUE} computes the trace of the smoothing matrix, the estimate of the error's variance and 
#'        the Generalized Cross Validation parameter, for value of \code{lambda}.
#' @param CPP_CODE if \code{TRUE} avoids the computation of some additional elements, not necessary if the 
#'        functions working with the FEM basis are called with the flag \code{CPP_CODE=TRUE}
#' @return A list with the following variables:
#'          \item{\code{FELSPLOBJ}}{A FEM object of the FEM type defined by the coefficients vector resulting from smoothing.}
#'          \item{\code{LAPLACEFD}}{A FEM object of the FEM type for the value of the Laplace operator}
#'          \item{\code{DOF}}{If GCV is \code{TRUE}, a vector with the trace of the smoothing matrix for each lambda.}
#'          \item{\code{sigma}}{If GCV is \code{TRUE}, a vector with the estimate of the standard deviation of the error for each lambda.}
#'          \item{\code{GCV}}{If GCV is \code{TRUE}, a vector with the GCV index for each lambda.}
#' @description Compute a solution for a Spatial Spline problem following the model in: Sangalli, Ramsay, Ramsay (2013).
#' @usage smooth.FEM.basis(locations = NULL, observations, basisobj, lambda, covariates = NULL, BC = NULL, GCV = TRUE, CPP_CODE = TRUE)
#' @examples library(FEMr)
#' data(MeuseData)
#' data(MeuseBorder)
#' order=1
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = order)
#' plot(mesh)
#' data = log(MeuseData[,7])
#' basisobj = create.FEM.basis(mesh, order)
#' lambda = 10^3.5
#' ZincMeuse = smooth.FEM.basis(observations = data, basisobj = basisobj, lambda = lambda)
#' plot(ZincMeuse$felsplobj)

smooth.FEM.basis<-function(locations = NULL, observations, basisobj, lambda, covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
{
  bigsol = NULL  
  lambda = as.vector(lambda)
  
  if(CPP_CODE == FALSE)
  {
    print('R Code Execution')
    if(!is.null(BC))
    {
      print('Dirichlet Boundary Conditions are ignored when CPP_CODE = FALSE. 
            If you want to use Dirichlet boundary conditions, please set CPP_CODE = TRUE')
    }
    
    bigsol = R_smooth.FEM.basis(locations, observations, basisobj, lambda, covariates, GCV)   
  }else
  {
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.basis(locations, observations, basisobj, lambda, covariates, BC, GCV)
  }
  
  numnodes = nrow(basisobj$mesh$nodes)
  
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  fit.FEM  = FEM(f, basisobj)
  PDEmisfit.FEM = FEM(g, basisobj)  
  
  reslist = NULL
  beta = getBetaCoefficients(locations, observations, fit.FEM, covariates)
  if(GCV == TRUE)
  {
    seq=getGCV(locations = locations, observations = observations, fit.FEM = fit.FEM, covariates = covariates, edfseq = bigsol[[2]])
    reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta, edf = bigsol[[2]], stderrseq = seq$stderrseq, GCVseq = seq$GCVseq)
  }else{
    reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta)
  }
  
  return(reslist)
}

#' Compute a solution for a Spatial Spline problem with PDE penalization
#' 
#' @param observations A vector specifying the observed values on the domain. 
#' The locations of the observations can be specified with the \code{locations} 
#' argument, otherwise the locations are intented to be the corresponding nodes of the mesh. 
#' \code{NA} values are admissible to indicate the missing value on the corresponding node.
#' @param locations A 2 column matrix where each row specifies the coordinates of the corresponding observation.
#' @param basisobj An an object of type FEM; See \code{\link{create.FEM.basis}}.
#' @param lambda A scalar smoothing parameter.
#' @param PDE_parameters A list containing the parameters of the penalizing PDE, with: \code{K} a 2-by-2 matrix indicating the diffusion coefficient matrix, \code{beta} a vector of length 2 with the coefficients of the advection coefficients and \code{c} a numeric indicating the reaction coefficient.
#' @param covariates A design matrix where each row represents the covariates associated to each row.
#' @param BC A list with two vectors: 
#'        \code{Indices}, a vector with the indices for the border points to apply a Dirichlet Border Condition;
#'        \code{Values} a vector with the values that the the nodes specified in \code{Indices} must assume.
#' @param GCV If \code{TRUE} computes the trace of the smoothing matrix, the estimate of the error's variance and 
#'        the Generalized Cross Validation parameter, for value of \code{lambda}.
#' @param CPP_CODE if \code{TRUE} avoids the computation of some additional elements, not necessary if the 
#'        functions working with the FEM basis are called with the flag \code{CPP_CODE=TRUE}
#' @return A list with the following variables:
#'          \item{\code{fit.FEM}}{A FEM object of the FEM type defined by the coefficients vector resulting from smoothing.}
#'          \item{\code{PDEmisfit.FEM}}{A FEM object of the FEM type for the value of the Laplace operator}
#'          \item{\code{edf}}{If GCV is \code{TRUE}, a vector with the trace of the smoothing matrix for each lambda.}
#'          \item{\code{sigma}}{If GCV is \code{TRUE}, a vector with the estimate of the standard deviation of the error for each lambda.}
#'          \item{\code{GCV}}{If GCV is \code{TRUE}, a vector with the GCV index for each lambda.}
#' @description Compute a solution for a for a Spatial Regression with PDE Penalization model.
#' @usage smooth.FEM.PDE.basis(locations = NULL, observations, basisobj, 
#'        lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = TRUE, 
#'        CPP_CODE = TRUE)
#' @examples 
#' library(FEMr)
#'
#' data(mesh.2D.simple)
#' plot(mesh)
#' observations = sin(pi*mesh$nodes[,1]) + rnorm(n = nrow(mesh$nodes), sd = 0.1)
#' 
#' basisobj = create.FEM.basis(mesh, 2)
#' 
#' # Smoothing coefficients
#' lambda = c(10^-2, 10^-1, 0.5, 5, 10)
#' 
#' # Anysotropic smoothing
#' PDE_parameters_anys = list(K = matrix(c(0.01,0,0,1), nrow = 2), b = c(0,0), c = 0)
#' FEM_CPP_PDE = smooth.FEM.PDE.basis(observations = observations, 
#'                                    basisobj = basisobj, lambda = lambda, PDE_parameters = PDE_parameters_anys)
#' plot(FEM_CPP_PDE$fit.FEM)
smooth.FEM.PDE.basis<-function(locations = NULL, observations, basisobj, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
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
  
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  fit.FEM  = FEM(f, basisobj)
  PDEmisfit.FEM = FEM(g, basisobj)  
  beta = getBetaCoefficients(locations, observations, fit.FEM, covariates)
  
  reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta=beta, edf = bigsol[[2]])
  return(reslist)
}


#' Compute a solution for a Spatial Spline problem with PDE penalization
#' 
#' @param observations A vector specifying the observed values on the domain. 
#' The locations of the observations can be specified with the \code{locations} 
#' argument, otherwise the locations are intented to be the corresponding nodes of the mesh. 
#' \code{NA} values are admissible to indicate the missing value on the corresponding node.
#' @param locations A 2 column matrix where each row specifies the coordinates of the corresponding observation.
#' @param basisobj An an object of type FEM; See \code{\link{create.FEM.basis}}.
#' @param lambda A scalar smoothing parameter.
#' @param PDE_parameters A list containing the space varying parameters of the penalizing PDE, with: \code{K} a function that for each point of the domain specified as a vector of length 2, returns  a 2-by-2 matrix indicating the diffusion coefficient matrix, \code{beta} a function that for each point of the domain specified as a vector of length 2, returns  a vector of length 2 indicating the advection coefficients and \code{c} a function that for each point of the domain specified as a vector of length 2, returns  a numeric indicating the reacttion coefficient.
#' @param covariates A design matrix where each row represents the covariates associated to each row.
#' @param BC A list with two vectors: 
#'        \code{Indices}, a vector with the indices for the border points to apply a Dirichlet Border Condition;
#'        \code{Values} a vector with the values that the the nodes specified in \code{Indices} must assume.
#' @param GCV If \code{TRUE} computes the trace of the smoothing matrix, the estimate of the error's variance and 
#'        the Generalized Cross Validation parameter, for value of \code{lambda}.
#' @param CPP_CODE if \code{TRUE} avoids the computation of some additional elements, not necessary if the 
#'        functions working with the FEM basis are called with the flag \code{CPP_CODE=TRUE}
#' @return A list with the following variables:
#'          \item{\code{fit.FEM}}{A FEM object of the FEM type defined by the coefficients vector resulting from smoothing.}
#'          \item{\code{PDEmisfit.FEM}}{A FEM object of the FEM type for the value of the Laplace operator}
#'          \item{\code{edf}}{If GCV is \code{TRUE}, a vector with the trace of the smoothing matrix for each lambda.}
#'          \item{\code{sigma}}{If GCV is \code{TRUE}, a vector with the estimate of the standard deviation of the error for each lambda.}
#'          \item{\code{GCV}}{If GCV is \code{TRUE}, a vector with the GCV index for each lambda.}
#' @description Compute a solution for a for a Spatial Regression with PDE Penalization model. The PDE's parameter are space-variant functions.
#' @usage smooth.FEM.PDE.SV.basis(locations = NULL, observations, basisobj, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = TRUE, CPP_CODE = TRUE)
#' @examples 
#' library(FEMr)
#' data(mesh.2D.rectangular)
#' observations = sin(0.2*pi*mesh$nodes[,1]) + rnorm(n = nrow(mesh$nodes), sd = 0.1)
#' basisobj = create.FEM.basis(mesh, 2)
#' # Smoothing coefficient
#' lambda = c(10^-2)
#' K_func<-function(points)
#' {
#' mat<-c(0.01,0,0,1)
#' as.vector(0.5*mat %*% t(points[,1]^2))
#' }
#' b_func<-function(points)
#' {
#' rep(c(0,0), nrow(points))
#' }
#' 
#' c_func<-function(points)
#' {
#' rep(c(0), nrow(points))
#' }
#' 
#' u_func<-function(points)
#' {
#' rep(c(0), nrow(points))
#' }
#' # Space-varying smoothing
#' PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)
#' FEM_CPP_PDE = smooth.FEM.PDE.SV.basis(observations = observations, 
#'              basisobj = basisobj, lambda = lambda, PDE_parameters = PDE_parameters)
#' plot(FEM_CPP_PDE$fit.FEM)
smooth.FEM.PDE.SV.basis<-function(locations = NULL, observations, basisobj, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
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
  
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  fit.FEM  = FEM(f, basisobj)
  PDEmisfit.FEM = FEM(g, basisobj)  
  beta = getBetaCoefficients(locations, observations, fit.FEM, covariates)
  
  reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta, edf = bigsol[[2]])
  return(reslist)
}



getBetaCoefficients<-function(locations, observations, fit.FEM, covariates)
{
  loc_nodes = NULL
  fnhat = NULL
  betahat = NULL
  
  if(!is.null(covariates))
  {
    if(is.null(locations))
    {
      loc_nodes = (1:length(observations))[!is.na(observations)]
      fnhat = fit.FEM$coefmat[loc_nodes,]
    }else{
      loc_nodes = 1:length(observations)
      fnhat = eval.FEM(FEM = fit.FEM, locations = locations, CPP_CODE = FALSE)
    }
    
    betahat = lm.fit(covariates,observations-fnhat)$coefficients
  }
  
 return(betahat)
}


getGCV<-function(locations, observations, fit.FEM, covariates = NULL, edfseq)
{
  loc_nodes = NULL
  fnhat = NULL
  
  edfseq = as.vector(edfseq)
    
  if(is.null(locations))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    fnhat = fit.FEM$coefmat[loc_nodes,]
  }else{
    loc_nodes = 1:length(observations)
    fnhat = eval.FEM(FEM = fit.FEM, locations = locations, CPP_CODE = FALSE)
  }
  
  zhat = NULL
  zhat = matrix(nrow = length(loc_nodes), ncol = length(edfseq))
  if(!is.null(covariates))
  {
    desmatprod = ( solve( t(covariates) %*% covariates ) ) %*% t(covariates)
    for ( i in 1:length(edfseq))
    {
      betahat  = desmatprod %*% (observations-fnhat)[,i]
      zhat[,i]     = covariates %*% betahat + fnhat[,i]
    }
  }else{
    zhat = fnhat
  }
  
  np = length(loc_nodes)
  
  stderr2seq = numeric(length(edfseq))
  GCVseq       = numeric(length(edfseq))
  
  zhat <- as.matrix(zhat)
  
  for (i in 1:length(edfseq))
  {
    stderr2seq[i] = t(observations[loc_nodes] - zhat[,i]) %*% (observations[loc_nodes] - zhat[,i]) / ( np - edfseq[i] )
    GCVseq[i] = ( np / ( np - edfseq[i] )) * stderr2seq[i]
  }
  
  return(list(stderrseq = sqrt(stderr2seq), GCVseq = GCVseq))
}
