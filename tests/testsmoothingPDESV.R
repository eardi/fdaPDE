## This script tests 
## - PDE smoothing with SV coefficients 
## - 1st order FEs
## - C++ code

library(FEMr)

data(mesh.2D.rectangular)
observations = sin(0.2*pi*mesh$nodes[,1]) + rnorm(n = nrow(mesh$nodes), sd = 0.1)

basisobj = create.FEM.basis(mesh, 2)

lambda = c(10^-2)

K_func<-function(points)
{
  mat<-c(0.01,0,0,1)
  as.vector(0.5*mat %*% t(points[,1]^2))
}

b_func<-function(points)
{
  rep(c(0,0), nrow(points))
}

c_func<-function(points)
{
  rep(c(0), nrow(points))
}
u_func<-function(points)
{
  rep(c(0), nrow(points))
}
# Space-varying smoothing
PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)
FEM_CPP_PDE_SV = smooth.FEM.PDE.SV.basis(observations = observations, 
                                      basisobj = basisobj, lambda = lambda, PDE_parameters = PDE_parameters)
print(FEM_CPP_PDE_SV$fit.FEM$coefmat)
