#' Compute a spatial regression model with differential regularization: stationary and isotropic case
#' 
#' @param observations A #observations vector with the observed values on the domain. 
#' The locations of the observations can be specified with the \code{locations} 
#' argument, otherwise the locations are intented to be the corresponding nodes of the mesh****. 
#' \code{NA} values are admissible to indicate that the node is not associated with any observed data value.
#' @param locations A #observations-by-2 matrix where each row specifies the spatial coordinates of the corresponding observation in \code{observations}. ***
#' @param basisobj A FEM object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values the spatial field must take at the nodes indicated in \code{BC_Indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, corresponding the values of the smoothing parameter specified in \code{lambda}.
#' @param CPP_CODE Boolean. If \code{TRUE} it avoids the computation of some additional quantities, that are not necessary if the 
#'        functions using the FEM basis are called with the flag \code{CPP_CODE=TRUE}
#' @return A list with the following variables:
#' \item{\code{fit.FEM}}{A FEM object that represents the fitted spatial field.}
#' \item{\code{PDEmisfit.FEM}}{A FEM object that represents the Laplacian of the estimated spatial field.}
#' \item{\code{beta}}{If covariates is not \code{NULL}, a #covariates vector with the regression coefficients associated with each covariate.}
#' \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#' \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#' \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each for each value of the smoothing parameter specified in \code{lambda}.}
#' @description This function implements a spatial regression model with differential regularization; isotropic and stationary case. In particular, the regularizing term involves the Laplacian of the spatial field. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
#' @usage smooth.FEM.basis(locations = NULL, observations, basisobj, lambda, 
#'        covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
#' @references Sangalli, L.M., Ramsay, J.O. & Ramsay, T.O., 2013. Spatial spline regression models. Journal of the Royal Statistical Society. Series B: Statistical Methodology, 75(4), pp. 681-703.
#' @examples
#' ## Upload the Meuse data and a domain boundary for these data
#' data(MeuseData)
#' data(MeuseBorder)
#' ## Create a triangular mesh for these data with the provided boundary
#' order=1
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = order)
#' plot(mesh)
#' data = log(MeuseData[,7])
#' basisobj = create.FEM.basis(mesh, order)
#' lambda = 10^3.5
#' ZincMeuse = smooth.FEM.basis(observations = data, basisobj = basisobj, lambda = lambda)
#' plot(ZincMeuse$fit.FEM)

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
    seq=getGCV(locations = locations, observations = observations, fit.FEM = fit.FEM, covariates = covariates, edf = bigsol[[2]])
    reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta, edf = bigsol[[2]], stderr = seq$stderr, GCV = seq$GCV)
  }else{
    reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta)
  }
  
  return(reslist)
}

#' Compute a spatial regression model with differential regularization: anysotropic case (elliptic PDE)
#' 
#' @param observations A #observations vector with the observed values on the domain. 
#' The locations of the observations can be specified with the \code{locations} 
#' argument, otherwise the locations are intented to be the corresponding nodes of the mesh****. 
#' \code{NA} values are admissible to indicate that the node is not associated with any observed data value.
#' @param locations A #observations-by-2 matrix where each row specifies the spatial coordinates of the corresponding observation in \code{observations}. ***
#' @param basisobj A FEM object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param PDE_parameters A list specifying the parameters of the elliptic PDE in the regularizing term: \code{K}, the 2-by-2 matrix of the diffusion tensor; \code{beta}, a 2entries vector of length 2 with the coefficients of the advection coefficients and \code{c} a numeric indicating the reaction coefficient.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values the spatial field must take at the nodes indicated in \code{BC_Indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, corresponding the values of the smoothing parameter specified in \code{lambda}.
#' @param CPP_CODE Boolean. If \code{TRUE} it avoids the computation of some additional quantities, that are not necessary if the 
#'        functions using the FEM basis are called with the flag \code{CPP_CODE=TRUE}
#' @return A list with the following variables:
#'          \item{\code{fit.FEM}}{A FEM object that represents the fitted spatial field.}
#'          \item{\code{PDEmisfit.FEM}}{A FEM object that represents the PDE misfit for the estimated spatial field.}
#'          \item{\code{beta}}{If covariates is not \code{NULL}, a #covariates vector with the regression coefficients associated with each covariate.}
#'          \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each for each value of the smoothing parameter specified in \code{lambda}.}
#' @description This function implements a spatial regression model with differential regularization; anysotropic case. In particular, the regularizing term involves a second order elliptic PDE, that governs the phenomenon behavior. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
#' @usage smooth.FEM.PDE.basis(locations = NULL, observations, basisobj, 
#'        lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, 
#'        CPP_CODE = TRUE)
#' @references Azzimonti, L., Sangalli, L.M., Secchi, P., Domanin, M., and Nobile, F., 2014. Blood flow velocity field estimation via spatial regression with PDE penalization Blood flow velocity field estimation via spatial regression with PDE penalization. DOI. 10.1080/01621459.2014.946036. Available at: http://dx.doi.org/10.1080/01621459.2014.946036
#'  Azzimonti, L., Nobile, F., Sangalli, L.M., and Secchi, P., 2014. Mixed Finite Elements for Spatial Regression with PDE Penalization. SIAM/ASA Journal on Uncertainty Quantification, 2(1), pp.305-335. Available at: http://epubs.siam.org/doi/abs/10.1137/130925426.
#' @examples 
#' data(mesh.2D.simple)
#' plot(mesh.2D.simple)
#' observations = sin(pi*mesh.2D.simple$nodes[,1]) + rnorm(n = nrow(mesh.2D.simple$nodes), sd = 0.1)
#' 
#' basisobj = create.FEM.basis(mesh.2D.simple, 2)
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
  
  reslist = NULL
  beta = getBetaCoefficients(locations, observations, fit.FEM, covariates)
  if(GCV == TRUE)
  {
    seq=getGCV(locations = locations, observations = observations, fit.FEM = fit.FEM, covariates = covariates, edf = bigsol[[2]])
    reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta, edf = bigsol[[2]], stderr = seq$stderr, GCV = seq$GCV)
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
#' @param basisobj An object of class FEM; See \code{\link{create.FEM.basis}}.
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
#'          \item{\code{fit.FEM}}{A FEM object of the FEM class defined by the coefficients vector resulting from smoothing.}
#'          \item{\code{PDEmisfit.FEM}}{A FEM object of the FEM class for the value of the Laplace operator}
#'          \item{\code{beta}}{If covariates is not \code{NULL}, a vector with the linear coefficients associated with each covariate.}
#'          \item{\code{edf}}{If GCV is \code{TRUE}, a vector with the trace of the smoothing matrix for each penalization parameter in the vector \code{lambda}.}
#'          \item{\code{stderr}}{If GCV is \code{TRUE}, a vector with the estimate of the standard deviation of the error for each penalization parameter in the vector \code{lambda}.}
#'          \item{\code{GCV}}{If GCV is \code{TRUE}, a vector with the GCV index for each penalization parameter in the vector \code{lambda}.}
#' @description Compute a solution for a for a Spatial Regression with PDE Penalization model. The PDE's parameter are space-variant functions.
#' @usage smooth.FEM.PDE.SV.basis(locations = NULL, observations, basisobj, 
#'  lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, 
#'  CPP_CODE = TRUE)
#' @references Azzimonti, L. et al., 2014. Blood flow velocity field estimation via spatial regression with PDE penalization Blood flow velocity field estimation via spatial regression with PDE penalization. , (September), pp.37-41.
#'  Azzimonti, L. et al., 2014. Mixed Finite Elements for Spatial Regression with PDE Penalization. SIAM/ASA Journal on Uncertainty Quantification, 2(1), pp.305-335. Available at: http://epubs.siam.org/doi/abs/10.1137/130925426.
#' @examples 
#' library(FEMr)
#' data(mesh.2D.rectangular)
#' observations = sin(0.2*pi*mesh.2D.rectangular$nodes[,1]) + 
#' rnorm(n = nrow(mesh.2D.rectangular$nodes), sd = 0.1)
#' basisobj = create.FEM.basis(mesh.2D.rectangular, 2)
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
  
  reslist = NULL
  beta = getBetaCoefficients(locations, observations, fit.FEM, covariates)
  if(GCV == TRUE)
  {
    seq=getGCV(locations = locations, observations = observations, fit.FEM = fit.FEM, covariates = covariates, edf = bigsol[[2]])
    reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta, edf = bigsol[[2]], stderr = seq$stderr, GCV = seq$GCV)
  }else{
    reslist=list(fit.FEM=fit.FEM,PDEmisfit.FEM=PDEmisfit.FEM, beta = beta)
  }
  
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


getGCV<-function(locations, observations, fit.FEM, covariates = NULL, edf)
{
  loc_nodes = NULL
  fnhat = NULL
  
  edf = as.vector(edf)
    
  if(is.null(locations))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    fnhat = fit.FEM$coefmat[loc_nodes,]
  }else{
    loc_nodes = 1:length(observations)
    fnhat = eval.FEM(FEM = fit.FEM, locations = locations, CPP_CODE = FALSE)
  }
  
  zhat = NULL
  zhat = matrix(nrow = length(loc_nodes), ncol = length(edf))
  if(!is.null(covariates))
  {
    desmatprod = ( solve( t(covariates) %*% covariates ) ) %*% t(covariates)
    for ( i in 1:length(edf))
    {
      betahat  = desmatprod %*% (observations-fnhat)[,i]
      zhat[,i]     = covariates %*% betahat + fnhat[,i]
    }
  }else{
    zhat = fnhat
  }
  
  np = length(loc_nodes)
  
  stderr2 = numeric(length(edf))
  GCV       = numeric(length(edf))
  
  zhat <- as.matrix(zhat)
  
  for (i in 1:length(edf))
  {
    stderr2[i] = t(observations[loc_nodes] - zhat[,i]) %*% (observations[loc_nodes] - zhat[,i]) / ( np - edf[i] )
    GCV[i] = ( np / ( np - edf[i] )) * stderr2[i]
  }
  
  return(list(stderr = sqrt(stderr2), GCV = GCV))
}
