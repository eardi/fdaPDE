library(FEMr)

data(mesh.2D.simple)
plot(mesh)
observations = sin(pi*mesh$nodes[,1]) + rnorm(n = nrow(mesh$nodes), sd = 0.1)

basisobj = create.FEM.basis(mesh, 2)

# Smoothing coefficients
lambda = c(10^-2, 10^-1, 0.5, 5, 10)

# Anysotropic smoothing
PDE_parameters_anys = list(K = matrix(c(0.01,0,0,1), nrow = 2), b = c(0,0), c = 0)
FEM_CPP_PDE = smooth.FEM.PDE.basis(observations = observations, 
                                   basisobj = basisobj, lambda = lambda, PDE_parameters = PDE_parameters_anys)
plot(FEM_CPP_PDE$fit.FEM)
