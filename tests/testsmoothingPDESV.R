library(RPDE)
setwd("~/workspace/RPDE/RScripts")

order = 1
mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                     segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)), order = order)

basisobj = create.FEM.basis(mesh, order)

#  smooth the data without covariates
lambda = c(1,2,3)

## data diviso in due
locations = rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0))
observations = c(1,2,1,2,1)
data = c(1,2,1,2,1)
covariates = cbind(c(1, 2, 3, 4, 5))
BC = NULL

K_func<-function(points) {rep(c(1,0,0,1), nrow(points))}
beta_func<-function(points){rep(c(0,0), nrow(points))}
c_func<-function(points){rep(c(0), nrow(points))}
u_func<-function(points){rep(c(0), nrow(points))}
PDE_parameters = list(K = K_func, beta = beta_func, c = c_func, u = u_func)
output_CPP_PDE_SV = smooth.FEM.PDE.SV.basis(locations  = as.matrix(locations), observations = data, 
                                        basisobj = basisobj, lambda = lambda, PDE_parameters, covariates = covariates)

print(output_CPP_PDE_SV$felsplobj$coefmat)
