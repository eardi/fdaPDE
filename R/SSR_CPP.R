#dyn.load("../Release/FEMr.so")

CPP_smooth.FEM.basis<-function(locations, observations, basisobj, lambda, covariates = NULL, BC = NULL, GCV)
{
  # SMOOTH.FEM.FD Compute a solution for a Spatial Spline problem 
  #
  #     Arguments:
  # FELSPLOBJ a FELspline object.
  # LAMBDA    a scalar smoothing parameter.
  # DATA      a n-by-1 matrix, representing the
  #           set of noisy observations of the surface values.
  # DESMAT    a n-by-p matrix representing the design matrix of the
  #           regression
  # BINDEX    a vector with the nodes indexes where to apply Dirichlet condition
  #            if NULL (NO Dirichlet condition)
  # BVALUES   a vector containing the values to apply to the nodes specified with
  #           BINDEX
  #
  #     Output:
  # FELSPLOBJ  ...  A FD object of the FEM type defined by the coefficient
  #                 vector resulting from smoothing
  # LAPLACEFD  ...  A FD object of the FEM type for the value of the 
  #                 Laplace operator if order == 2, or empty if order == 1
  
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  basisobj$mesh$triangles = basisobj$mesh$triangles - 1
  basisobj$mesh$edges = basisobj$mesh$edges - 1
  basisobj$mesh$neighbors = basisobj$mesh$neighbors - 1
  
  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = 2)
  }
  
  if(is.null(BC$Indices))
  {
    BC$Indices<-vector(length=0)
  }else
  {
    BC$Indices<-as.vector(BC$Indices)-1
  }
  
  if(is.null(BC$Values))
  {
    BC$Values<-vector(length=0)
  }else
  {
    BC$Values<-as.vector(BC$Values)
  }
  
  ## Set propr type for correct C++ reading
  storage.mode(locations) <- "double"
  storage.mode(basisobj$mesh$points) <- "double"
  storage.mode(basisobj$mesh$triangles) <- "integer"
  storage.mode(basisobj$mesh$edges) <- "integer"
  storage.mode(basisobj$mesh$neighbors) <- "integer"
  storage.mode(basisobj$order) <- "integer"
  storage.mode(covariates) <- "double"
  storage.mode(lambda)<- "double"
  storage.mode(BC$Indices)<- "integer"
  storage.mode(BC$Values)<-"double"
  
  GCV = as.integer(GCV)
  storage.mode(GCV)<-"integer"
  
  ## Call C++ function
  bigsol <- .Call("regression_Laplace", locations, observations, basisobj$mesh, 
                  basisobj$order, lambda, covariates,
                  BC$Indices, BC$Values, GCV,
                  package = "FEMr")
  
  ## Reset them correctly
  #fdobj$basis$params$mesh$triangles = fdobj$basis$params$mesh$triangles + 1
  #fdobj$basis$params$mesh$edges = fdobj$basis$params$mesh$edges + 1
  #fdobj$basis$params$mesh$neighbors = fdobj$basis$params$mesh$neighbors + 1
  return(bigsol)
}

CPP_smooth.FEM.PDE.basis<-function(locations, observations, basisobj, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV)
{
  # SMOOTH.FEM.FD Compute a solution for a Spatial Spline problem 
  #
  #     Arguments:
  # FELSPLOBJ a FELspline object.
  # LAMBDA    a scalar smoothing parameter.
  # DATA      a n-by-1 matrix, representing the
  #           set of noisy observations of the surface values.
  # DESMAT    a n-by-p matrix representing the design matrix of the
  #           regression
  # BINDEX    a vector with the nodes indexes where to apply Dirichlet condition
  #            if NULL (NO Dirichlet condition)
  # BVALUES   a vector containing the values to apply to the nodes specified with
  #           BINDEX
  #
  #     Output:
  # FELSPLOBJ  ...  A FD object of the FEM type defined by the coefficient
  #                 vector resulting from smoothing
  # LAPLACEFD  ...  A FD object of the FEM type for the value of the 
  #                 Laplace operator if order == 2, or empty if order == 1
  
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  basisobj$mesh$triangles = basisobj$mesh$triangles - 1
  basisobj$mesh$edges = basisobj$mesh$edges - 1
  basisobj$mesh$neighbors = basisobj$mesh$neighbors - 1
  
  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = 2)
  }
  
  if(is.null(BC$Indices))
  {
    BC$Indices<-vector(length=0)
  }else
  {
    BC$Indices<-as.vector(BC$Indices)-1
  }
  
  if(is.null(BC$Values))
  {
    BC$Values<-vector(length=0)
  }else
  {
    BC$Values<-as.vector(BC$Values)
  }
  
  ## Set propr type for correct C++ reading
  storage.mode(locations) <- "double"
  storage.mode(basisobj$mesh$points) <- "double"
  storage.mode(basisobj$mesh$triangles) <- "integer"
  storage.mode(basisobj$mesh$edges) <- "integer"
  storage.mode(basisobj$mesh$neighbors) <- "integer"
  storage.mode(basisobj$order) <- "integer"
  storage.mode(covariates) <- "double"
  storage.mode(lambda)<- "double"
  storage.mode(BC$Indices)<- "integer"
  storage.mode(BC$Values)<-"double"
  storage.mode(GCV)<-"integer"
  
  storage.mode(PDE_parameters$K)<-"double"
  storage.mode(PDE_parameters$beta)<-"double"
  storage.mode(PDE_parameters$c)<-"double"
  
  ## Call C++ function
  bigsol <- .Call("regression_PDE", locations, observations, basisobj$mesh, 
                  basisobj$order, lambda, PDE_parameters$K, PDE_parameters$beta, PDE_parameters$c, covariates,
                  BC$Indices, BC$Values, GCV,
                  package = "FEMr")
  
  ## Reset them correctly
  #fdobj$basis$params$mesh$triangles = fdobj$basis$params$mesh$triangles + 1
  #fdobj$basis$params$mesh$edges = fdobj$basis$params$mesh$edges + 1
  #fdobj$basis$params$mesh$neighbors = fdobj$basis$params$mesh$neighbors + 1
  return(bigsol)
}

CPP_smooth.FEM.PDE.SV.basis<-function(locations, observations, basisobj, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV)
{
  
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  basisobj$mesh$triangles = basisobj$mesh$triangles - 1
  basisobj$mesh$edges = basisobj$mesh$edges - 1
  basisobj$mesh$neighbors = basisobj$mesh$neighbors - 1
  
  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = 2)
  }
  
  if(is.null(BC$Indices))
  {
    BC$Indices<-vector(length=0)
  }else
  {
    BC$Indices<-as.vector(BC$Indices)-1
  }
  
  if(is.null(BC$Values))
  {
    BC$Values<-vector(length=0)
  }else
  {
    BC$Values<-as.vector(BC$Values)
  }
  
  
  PDE_param_eval = NULL
  points_eval = matrix(CPP_get_evaluations_points(mesh = basisobj$mesh, order = basisobj$order),ncol = 2);
  PDE_param_eval$K = (PDE_parameters$K)(points_eval)
  PDE_param_eval$beta = (PDE_parameters$beta)(points_eval)
  PDE_param_eval$c = (PDE_parameters$c)(points_eval)
  PDE_param_eval$u = (PDE_parameters$u)(points_eval)
  
  #print(PDE_param_eval)
  ## Set propr type for correct C++ reading
  storage.mode(locations) <- "double"
  storage.mode(basisobj$mesh$points) <- "double"
  storage.mode(basisobj$mesh$triangles) <- "integer"
  storage.mode(basisobj$mesh$edges) <- "integer"
  storage.mode(basisobj$mesh$neighbors) <- "integer"
  storage.mode(basisobj$order) <- "integer"
  storage.mode(covariates) <- "double"
  storage.mode(lambda)<- "double"
  storage.mode(BC$Indices)<- "integer"
  storage.mode(BC$Values)<-"double"
  storage.mode(GCV)<-"integer"
  
  storage.mode(PDE_param_eval$K)<-"double"
  storage.mode(PDE_param_eval$beta)<-"double"
  storage.mode(PDE_param_eval$c)<-"double"
  storage.mode(PDE_param_eval$u)<-"double"
  
  ## Call C++ function
  bigsol <- .Call("regression_PDE_space_varying", locations, observations, basisobj$mesh, 
                  basisobj$order, lambda, PDE_param_eval$K, PDE_param_eval$beta, PDE_param_eval$c, PDE_param_eval$u, covariates,
                  BC$Indices, BC$Values, GCV,
                  package = "FEMr")
  
  ## Reset them correctly
  #fdobj$basis$params$mesh$triangles = fdobj$basis$params$mesh$triangles + 1
  #fdobj$basis$params$mesh$edges = fdobj$basis$params$mesh$edges + 1
  #fdobj$basis$params$mesh$neighbors = fdobj$basis$params$mesh$neighbors + 1
  return(bigsol)
}

CPP_eval.FEM.fd = function(X,Y,fdobj,fast)
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
  
  ## C++ indexes starts from zero, needed conversion.
  fdobj$basis$params$mesh$triangles = fdobj$basis$params$mesh$triangles - 1
  fdobj$basis$params$mesh$edges = fdobj$basis$params$mesh$edges - 1
  fdobj$basis$params$mesh$neighbors = fdobj$basis$params$mesh$neighbors - 1
  
  # Imposing types, this is necessary for correct reading from C++
  storage.mode(fdobj$basis$params$mesh$points) <- "double"
  storage.mode(fdobj$basis$params$mesh$triangles) <- "integer"
  storage.mode(fdobj$basis$params$mesh$edges) <- "integer"
  storage.mode(fdobj$basis$params$mesh$neighbors) <- "integer"
  storage.mode(fdobj$basis$params$order) <- "integer"
  fdobj$coef = as.vector(fdobj$coef)
  storage.mode(fdobj$coef) <- "double"
  storage.mode(X) <- "double"
  storage.mode(Y) <- "double"
  storage.mode(fast) <- "integer"
  
  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  evalmat <- .Call("eval_FEM_fd",fdobj$basis$params$mesh, as.vector(X), as.vector(Y), as.vector(fdobj$coef), 
                   fdobj$basis$params$order, fast,
                   package = "FEMr")
  
  #Returning the evaluation matrix
  evalmat
}

CPP_get_evaluations_points = function(mesh, order)
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
  
  ## C++ indexes starts from zero, needed conversion.
  mesh$triangles = mesh$triangles - 1
  mesh$edges = mesh$edges - 1
  mesh$neighbors = mesh$neighbors - 1
  
  # Imposing types, this is necessary for correct reading from C++
  storage.mode(mesh$points) <- "double"
  storage.mode(mesh$triangles) <- "integer"
  storage.mode(mesh$edges) <- "integer"
  storage.mode(mesh$neighbors) <- "integer"
  storage.mode(order) <- "integer"
  
  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  points <- .Call("get_integration_points",mesh, order,
                  package = "FEMr")
  
  #Returning the evaluation matrix
  points
}