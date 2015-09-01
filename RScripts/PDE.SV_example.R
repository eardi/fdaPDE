library(FEMr)
setwd("~/git/FEMr/RScripts")
# nodes = rbind(c(1,1), c(20,1), c(20,5), c(1,5),cbind(runif(20, min = 1, max = 20), runif(20, min = 1, max = 5)))
# mesh = create.MESH.2D(nodes, segments = rbind(c(1,2),c(2,3),c(3,4),c(4,1)))
# mesh = refine.MESH.2D(mesh, maximum_area = 0.1, minimum_angle = 30)
# plot(mesh)

#save(mesh,file = "mesh.2D.long")
#data(file = "mesh.2D.long")

data(mesh.2D.rectangular)
observations = sin(0.2*pi*mesh$nodes[,1]) + rnorm(n = nrow(mesh$nodes), sd = 0.1)

basisobj = create.FEM.basis(mesh, 2)


# Smoothing coefficients
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
FEM_CPP_PDE = smooth.FEM.PDE.SV.basis(observations = observations, 
                                   basisobj = basisobj, lambda = lambda, PDE_parameters = PDE_parameters)
plot(FEM_CPP_PDE$fit.FEM)

