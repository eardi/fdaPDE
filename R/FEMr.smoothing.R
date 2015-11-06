#' Spatial regression with differential regularization: stationary and isotropic case (Laplacian)
#' 
#' @param observations A #observations vector with the observed data values over the domain. 
#' The locations of the observations can be specified with the \code{locations} 
#' argument, otherwise the locations are intented to be the corresponding nodes of the mesh****. 
#' \code{NA} values are admissible to indicate that the node is not associated with any observed data value.
#' @param locations A #observations-by-2 matrix where each row specifies the spatial coordinates of the corresponding observation in \code{observations}. ***
#' @param FEMbasis A FEM object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_Indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @param CPP_CODE Boolean. If \code{TRUE} it avoids the computation of some additional quantities, that are not necessary if the 
#'        functions using the FEM basis are called with the flag \code{CPP_CODE=TRUE}.
#' @return A list with the following variables:
#' \item{\code{fit.FEM}}{A FEM object that represents the fitted spatial field.}
#' \item{\code{PDEmisfit.FEM}}{A FEM object that represents the Laplacian of the estimated spatial field.}
#' \item{\code{beta}}{If covariates is not \code{NULL}, a vector of length #covariates with the regression coefficients associated with each covariate.}
#' \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#' \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#' \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' @description This function implements a spatial regression model with differential regularization; isotropic and stationary case. In particular, the regularizing term involves the Laplacian of the spatial field. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
#' @usage smooth.FEM.basis(locations = NULL, observations, FEMbasis, lambda, 
#'        covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
#' @references Sangalli, L.M., Ramsay, J.O. & Ramsay, T.O., 2013. Spatial spline regression models. Journal of the Royal Statistical Society. Series B: Statistical Methodology, 75(4), pp. 681-703.
#' @examples
#' library(FEMr)
#' ## Upload the Meuse data and a domain boundary for these data
#' data(MeuseData)
#' data(MeuseBorder)
#' ## Create a triangular mesh for these data with the provided boundary and plot it
#' order=1
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = order)
#' plot(mesh)
#' ## Create the Finite Element basis 
#' FEMbasis = create.FEM.basis(mesh, order)
#' ## Estimate ***** using as covariates ****, setting the smoothing parameter to **** 10^3.5
#' data = log(MeuseData[,"zinc"])
#' lambda = 10^3.5
#' ZincMeuse = smooth.FEM.basis(observations = data, FEMbasis = FEMbasis, lambda = lambda)
#' ## Plot the estimated spatial field 
#' plot(ZincMeuse$fit.FEM)
#' # Now repeat the analysis using as covariates the square root of the log-distance 
#' # from river \code{sqrt(dist.log(m))} and the altitude \code{elev}
#' desmat = matrix(1,nrow=nrow(MeuseData),ncol=2)
#' desmat[,1] = sqrt(MeuseData[,"dist.log(m)"])
#' desmat[,2] = MeuseData[,"elev"]
#' ZincMeuseCovar = smooth.FEM.basis(observations = data, covariates = desmat, FEMbasis = FEMbasis, lambda = lambda)
#' # Plot of the non parametric part (f) of the regression model y_i = beta_1 x_i1 + beta_2 x_i2 + f
#' plot(ZincMeuseCovar$fit.FEM)
#' # Print covariates' regression coefficients
#' print(ZincMeuseCovar$beta)


smooth.FEM.basis<-function(locations = NULL, observations, FEMbasis, lambda, covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
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
    
    bigsol = R_smooth.FEM.basis(locations, observations, FEMbasis, lambda, covariates, GCV)   
  }else
  {
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.basis(locations, observations, FEMbasis, lambda, covariates, BC, GCV)
  }
  
  numnodes = nrow(FEMbasis$mesh$nodes)
  
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  fit.FEM  = FEM(f, FEMbasis)
  PDEmisfit.FEM = FEM(g, FEMbasis)  
  
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

#' Spatial regression with differential regularization: anysotropic case (elliptic PDE)
#' 
#' @param observations A #observations vector with the observed data values over the domain. 
#' The locations of the observations can be specified with the \code{locations} 
#' argument, otherwise the locations are intented to be the corresponding nodes of the mesh****. 
#' \code{NA} values are admissible to indicate that the node is not associated with any observed data value.
#' @param locations A #observations-by-2 matrix where each row specifies the spatial coordinates of the corresponding observation in \code{observations}. ***
#' @param FEMbasis A FEM object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param PDE_parameters A list specifying the parameters of the elliptic PDE in the regularizing term: \code{K}, a 2-by-2 matrix of diffusion coefficients; \code{beta}, a vector of length 2 of advection coefficients;  \code{c}, a scalar reaction coefficient.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_Indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @param CPP_CODE Boolean. If \code{TRUE} it avoids the computation of some additional quantities, that are not necessary if the 
#'        functions using the FEM basis are called with the flag \code{CPP_CODE=TRUE}.
#' @return A list with the following variables:
#'          \item{\code{fit.FEM}}{A FEM object that represents the fitted spatial field.}
#'          \item{\code{PDEmisfit.FEM}}{A FEM object that represents the PDE misfit for the estimated spatial field.}
#'          \item{\code{beta}}{If covariates is not \code{NULL}, a vector of length #covariates with the regression coefficients associated with each covariate.}
#'          \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' @description This function implements a spatial regression model with differential regularization; anysotropic case. In particular, the regularizing term involves a second order elliptic PDE, that models the space-variation of the phenomenon. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
#' @usage smooth.FEM.PDE.basis(locations = NULL, observations, FEMbasis, 
#'        lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, 
#'        CPP_CODE = TRUE)
#' @references Azzimonti, L., Sangalli, L.M., Secchi, P., Domanin, M., and Nobile, F., 2014. Blood flow velocity field estimation via spatial regression with PDE penalization Blood flow velocity field estimation via spatial regression with PDE penalization. DOI. 10.1080/01621459.2014.946036. 
#'  Azzimonti, L., Nobile, F., Sangalli, L.M., and Secchi, P., 2014. Mixed Finite Elements for Spatial Regression with PDE Penalization. SIAM/ASA Journal on Uncertainty Quantification, 2(1), pp.305-335. 
#' @examples 
#' data(mesh.2D.simple)
#' plot(mesh.2D.simple)
#' observations = sin(pi*mesh.2D.simple$nodes[,1]) + rnorm(n = nrow(mesh.2D.simple$nodes), sd = 0.1)
#' 
#' FEMbasis = create.FEM.basis(mesh.2D.simple, 2)
#' 
#' # Smoothing coefficients
#' lambda = c(10^-2, 10^-1, 0.5, 5, 10)
#' 
#' # Anysotropic smoothing
#' PDE_parameters_anys = list(K = matrix(c(0.01,0,0,1), nrow = 2), b = c(0,0), c = 0)
#' FEM_CPP_PDE = smooth.FEM.PDE.basis(observations = observations, 
#'                                    FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters_anys)
#' plot(FEM_CPP_PDE$fit.FEM)
smooth.FEM.PDE.basis<-function(locations = NULL, observations, FEMbasis, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
{
  bigsol = NULL  
  lambda = as.vector(lambda)
  
  if(CPP_CODE == FALSE)
  {
    print('Function implemented only in C++, turn CPP_CODE = TRUE')  
  }else
  {
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.basis(locations, observations, FEMbasis, lambda, PDE_parameters, covariates, BC, GCV)
  }
  
  numnodes = nrow(FEMbasis$mesh$nodes)
  
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  fit.FEM  = FEM(f, FEMbasis)
  PDEmisfit.FEM = FEM(g, FEMbasis)  
  
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


#' Spatial regression with differential regularization: anysotropic and non-stationary case (elliptic PDE with space-varying coefficients)
#' 
#' @param observations A #observations vector with the observed data values over the domain. 
#' The locations of the observations can be specified with the \code{locations} 
#' argument, otherwise the locations are intented to be the corresponding nodes of the mesh****. 
#' \code{NA} values are admissible to indicate that the node is not associated with any observed data value.
#' @param locations A #observations-by-2 matrix where each row specifies the spatial coordinates of the corresponding observation in \code{observations}. ***
#' @param FEMbasis A FEM object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param PDE_parameters A list specifying the space-varying parameters of the elliptic PDE in the regularizing term: \code{K}, a function that for each spatial location in the spatial domain 
#' (indicated by the vector of the 2 spatial coordinates) returns a 2-by-2 matrix of diffusion coefficients; \code{beta}, a function that for each spatial location in the spatial domain returns a vector of length 2 of transport coefficients;  \code{c}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_Indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @param CPP_CODE Boolean. If \code{TRUE} it avoids the computation of some additional quantities, that are not necessary if the 
#'        functions using the FEM basis are called with the flag \code{CPP_CODE=TRUE}.
#' @return A list with the following variables:
#'          \item{\code{fit.FEM}}{A FEM object that represents the fitted spatial field.}
#'          \item{\code{PDEmisfit.FEM}}{A FEM object that represents the PDE misfit for the estimated spatial field.}
#'          \item{\code{beta}}{If covariates is not \code{NULL}, a vector of lenght #covariates with the regression coefficients associated with each covariate.}
#'          \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' @description This function implements a spatial regression model with differential regularization; anysotropic case. In particular, the regularizing term involves a second order elliptic PDE with space-varying coefficients, that models the space-variation of the phenomenon. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
#' @usage smooth.FEM.PDE.SV.basis(locations = NULL, observations, FEMbasis, 
#'  lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, 
#'  CPP_CODE = TRUE)
#' @references Azzimonti, L., Sangalli, L.M., Secchi, P., Domanin, M., and Nobile, F., 2014. Blood flow velocity field estimation via spatial regression with PDE penalization Blood flow velocity field estimation via spatial regression with PDE penalization. DOI. 10.1080/01621459.2014.946036. 
#'  Azzimonti, L., Nobile, F., Sangalli, L.M., and Secchi, P., 2014. Mixed Finite Elements for Spatial Regression with PDE Penalization. SIAM/ASA Journal on Uncertainty Quantification, 2(1), pp.305-335. 
#' @examples 
#' data(mesh.2D.rectangular)
#' FEMbasis = create.FEM.basis(mesh.2D.rectangular, 2)
#' observations = sin(0.2*pi*mesh.2D.rectangular$nodes[,1]) + 
#' rnorm(n = nrow(mesh.2D.rectangular$nodes), sd = 0.1)
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
#'              FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters)
#' plot(FEM_CPP_PDE$fit.FEM)
smooth.FEM.PDE.SV.basis<-function(locations = NULL, observations, FEMbasis, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
{
  bigsol = NULL  
  lambda = as.vector(lambda)
  
  if(CPP_CODE == FALSE)
  {
    print('Function implemented only in C++, turn CPP_CODE = TRUE')  
  }else
  {
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.SV.basis(locations, observations, FEMbasis, lambda, PDE_parameters, covariates, BC, GCV)
  }
  
  numnodes = nrow(FEMbasis$mesh$nodes)
  
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  fit.FEM  = FEM(f, FEMbasis)
  PDEmisfit.FEM = FEM(g, FEMbasis)  
  
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
      fnhat = fit.FEM$coeff[loc_nodes,]
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
    fnhat = fit.FEM$coeff[loc_nodes,]
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
