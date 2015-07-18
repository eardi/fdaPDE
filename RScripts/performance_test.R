library(rFEM)
#setwd("~/workspace/RPDE/RScripts")

order = 1
num_subdivisions = 40

x = seq(from = -100, to = 100, length.out = num_subdivisions)
y = seq(from = -100, to = 100, length.out = num_subdivisions)
nodes = expand.grid(x, y)
mesh<-create.MESH.2D(nodes=nodes, order = order)
#plot(mesh)

#mesh<-refine.MESH.2D(mesh,maximum_area = 0.05, delaunay = T)
#plot(mesh)

basisobj = create.FEM.basis(mesh, order)

lambda = 1

observation <-function(nodes)
{
  3*cos(0.2*nodes[,1])*cos(0.2*nodes[,2]) + rnorm(nrow(nodes), mean = 0, sd = 1)
}

data = observation(nodes)
BC = NULL

#BorderIndices = (1:length(mesh$nodesmarker))[mesh$nodesmarker==1]

#Indices= BorderIndices
#Values = rep(x = 0, times = length(Indices))
#BC = list(Indices = Indices, Values = Values)

## Due tipologie input
# output = smooth.LAPLACE.basis(locations  = locations, observations = data, 
#                                    basisobj = basisobj, lambda = lambda, covariates = covariates,
#                                    CPP_CODE = FALSE)

# output = smooth.LAPLACE.basis(locations  = as.matrix(locations), observations = data, 
#                               basisobj = basisobj, lambda = lambda, covariates = covariates,
#                               CPP_CODE = FALSE)


#OBS space varying smoothing function
K_func<-function(points)
{
  mat<-c(1,0,0,1)/300
  as.vector(mat %*% t(points[,1]^2+points[,2]^2))
}

beta_func<-function(points)
{
  rep(c(0,0), nrow(points))
}

c_func<-function(points)
{
  rep(c(0), nrow(points))
}
u_func<-function(points)
{
  rep(c(1), nrow(points))
}
PDE_parameters = list(K = K_func, beta = beta_func, c = c_func, u = u_func)
output = smooth.FEM.PDE.SV.basis(observations = data, 
                                        basisobj = basisobj, lambda = lambda, PDE_parameters = PDE_parameters, BC = BC)