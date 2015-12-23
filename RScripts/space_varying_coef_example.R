library(fdaPDE)
#setwd("~/workspace/RPDE/RScripts")

# Create a gid of points
x = seq(from = -100, to = 100, length.out = num_subdivisions)
y = seq(from = -100, to = 100, length.out = num_subdivisions)
nodes = expand.grid(x, y)

# Choose the order of the FE (either 1 or 2)
order = 2

# Create a mesh from the grid of points and plot it
mesh<-create.MESH.2D(nodes=nodes, order = order)
plot(mesh)

# If you want to refine the mesh, uncomment following two lines 
# (For very fine mesh consider working with order 1 FE , since the plotting is more efficient for FE of order 1)
# mesh<-refine.MESH.2D(mesh,maximum_area = 10, delaunay = T)
# plot(mesh)

FEMbasis = create.FEM.basis(mesh, order)

lambda = 1

observation <-function(nodes)
{
  3*cos(0.2*nodes[,1])*cos(0.2*nodes[,2]) + rnorm(nrow(nodes), mean = 0, sd = 1)
}

data = observation(nodes)

#create space varying smoothing coefficients (decading when far from zero)
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
# Group all coefficients in one object
PDE_parameters = list(K = K_func, beta = beta_func, c = c_func, u = u_func)

#Compute and print spatial field
FEM_CPP_PDE = smooth.FEM.PDE.sv.basis(observations = data,
                                      FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters)
plot(FEM_CPP_PDE$fit.FEM)
