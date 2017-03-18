rm(list = ls())
library(devtools)
library(ggplot2)
library(scales)
#library(fdaPDE)
load_all("/home/el425/store/git/fdaPDE/")
#load_all("/run/user/1000/gvfs/sftp:host=ssh.maths.cam.ac.uk,user=el425/store/DPMMS/el425/git/fdaPDE")
# settings
order = 1
lambda = 1
locations = NULL
covariates = NULL
BC = NULL

f = function(nodes) cos((nodes[,1]-0.5)*(pi))*cos((nodes[,2]-0.5)*(pi)) + rnorm(nrow(nodes),sd = 0.3)
PDE_parameters_anys = list(K = matrix(c(0.1,0,0,1), nrow = 2), b = c(0.1,0.1), c = 0)
K_func<-function(points)
{
  mat<-c(0.01,0,0,1)
  output = array(0, c(2, 2, nrow(points)))
  for (i in 1:nrow(points))
    output[,,i] = 0.5*mat %*% t(points[i,1]^2)
  output
}
b_func<-function(points)
{
  output = array(0, c(2, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = 0
  output
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
PDE_parameters_sv = list(K = K_func, b = b_func, c = c_func, u = u_func)

# x <- seq(0, 1, length.out = 4)
# y <- seq(0, 1, length.out = 4)
# nodes <- expand.grid(x = x, y = y)
# mesh<-create.MESH.2D(nodes=nodes, order = 1)
# FEMbasis = create.FEM.basis(mesh)
# observations = f(nodes)
# PDE_param = list(K = matrix(c(1,0,0,1), nrow = 2), b = c(0,0), c = 0)
# stiff_CPP = CPP_get.FEM.PDE.Matrix(FEMbasis, PDE_param)
# stiff_R = R_stiff(FEMbasis)
# PDE_param = list(K = matrix(c(0,0,0,0), nrow = 2), b = c(0,0), c = 1)
# mass_CPP = CPP_get.FEM.PDE.Matrix(FEMbasis, PDE_param)
# mass_R = R_mass(FEMbasis)

x <- seq(0, 1, length.out = 100)
y <- seq(0, 1, length.out = 100)
nodes <- expand.grid(x = x, y = y)
mesh<-create.MESH.2D(nodes=nodes, order = 1)
FEMbasis = create.FEM.basis(mesh)
observations = f(nodes)
oper_CPP = CPP_get.FEM.PDE.Matrix(FEMbasis, PDE_parameters_anys)

system.time(R_stiff(FEMbasis))
system.time(CPP_get.FEM.Mass.Matrix(FEMbasis))
system.time(CPP_get.FEM.PDE.Matrix(FEMbasis, PDE_parameters_anys))
system.time(CPP_get.FEM.PDE.sv.Matrix(FEMbasis, PDE_parameters_sv))
