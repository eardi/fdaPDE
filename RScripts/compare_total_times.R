rm(list = ls())
library(devtools)
library(ggplot2)
library(scales)
library(fdaPDE)
#load_all("/home/el425/store/git/fdaPDE/")
#load_all("/run/user/1000/gvfs/sftp:host=ssh.maths.cam.ac.uk,user=el425/store/DPMMS/el425/git/fdaPDE")
# settings
order = 1
lambda = 1
locations = NULL
covariates = NULL
BC = NULL

f = function(nodes) cos((nodes[,1]-0.5)*(pi))*cos((nodes[,2]-0.5)*(pi)) + rnorm(nrow(nodes),sd = 0.3)
PDE_parameters_anys = list(K = matrix(c(0.01,0,0,1), nrow = 2), b = c(0,0), c = 0)
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
# Simulation
times = NULL
for (sqrt_nnodes in seq(from = 10, by = 20, length.out = 7))
{
  for(rep in 1:10)
  {
    x <- seq(0, 1, length.out = sqrt_nnodes)
    y <- seq(0, 1, length.out = sqrt_nnodes)
    nodes <- expand.grid(x = x, y = y)
    mesh<-create.MESH.2D(nodes=nodes, order = order)
    FEMbasis = create.FEM.basis(mesh)
    observations = f(nodes)
    
    time_exec = system.time(smooth.FEM.basis(observations = observations, 
                                             FEMbasis = FEMbasis, lambda = lambda, 
                                             GCV = FALSE,
                                             CPP_CODE = FALSE))
    times = rbind(times,data.frame(rep = rep, sqrt_nnodes = sqrt_nnodes, nnodes = sqrt_nnodes^2, t_user = time_exec[1], t_system = time_exec[2], t_elapsed = time_exec[3], Code = "R Laplace"))
    
    time_exec = system.time(smooth.FEM.basis(observations = observations, 
                                             FEMbasis = FEMbasis, lambda = lambda, 
                                             GCV = FALSE,
                                             CPP_CODE = TRUE))
    times = rbind(times,data.frame(rep = rep, sqrt_nnodes = sqrt_nnodes, nnodes = sqrt_nnodes^2, t_user = time_exec[1], t_system = time_exec[2], t_elapsed = time_exec[3], Code = "C++ Laplace"))
    
    
    time_exec = system.time(smooth.FEM.PDE.basis(observations = observations, 
                                                 FEMbasis = FEMbasis, lambda = lambda, PDE_parameters =  PDE_parameters_anys,
                                                 GCV = FALSE,
                                                 CPP_CODE = TRUE))
    times = rbind(times,data.frame(rep = rep, sqrt_nnodes = sqrt_nnodes, nnodes = sqrt_nnodes^2, t_user = time_exec[1], t_system = time_exec[2], t_elapsed = time_exec[3], Code = "C++ PDE"))
    
    time_exec = system.time(smooth.FEM.PDE.sv.basis(observations = observations, 
                                                    FEMbasis = FEMbasis, lambda = lambda, PDE_parameters =  PDE_parameters_sv,
                                                    GCV = FALSE,
                                                    CPP_CODE = TRUE))
    times = rbind(times,data.frame(rep = rep, sqrt_nnodes = sqrt_nnodes, nnodes = sqrt_nnodes^2, t_user = time_exec[1], t_system = time_exec[2], t_elapsed = time_exec[3], Code = "C++ PDE-SV"))
  }
}

times_aggr = aggregate(. ~ nnodes+Code, data = times[,-1], FUN = function(x) c(mn = mean(x)))

plot <- ggplot(times_aggr, aes(nnodes, t_user, colour = Code)) + geom_point() +  geom_line() + ggtitle("Total execution times") + ylab("Execution time") + xlab("Number of nodes")
# log10 with exponents on tick labels
plot <- plot + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))
plot <- plot + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))

plot


x <- seq(0, 1, length.out = 2)
y <- seq(0, 1, length.out = 2)
nodes <- expand.grid(x = x, y = y)
mesh<-create.MESH.2D(nodes=nodes, order = 1)
FEMbasis = create.FEM.basis(mesh)
observations = f(nodes)
PDE_param = list(K = matrix(c(1,0,0,1), nrow = 2), b = c(0,0), c = 0)
stiff_CPP = CPP_get.FEM.PDE.Matrix(FEMbasis, PDE_param)
stiff_R = R_stiff(FEMbasis)

out = smooth.FEM.basis(observations = observations, 
                 FEMbasis = FEMbasis, lambda = lambda, 
                 GCV = FALSE,
                 CPP_CODE = TRUE)
plot(out$fit.FEM)

 # Rprof ( tf <- "log.log",  memory.profiling = TRUE )
 # output = smooth.FEM.basis(observations = observations,
 #                  FEMbasis = FEMbasis, lambda = lambda,
 #                  GCV = FALSE,
 #                  CPP_CODE = FALSE)
 # Rprof ( NULL ) ; print ( summaryRprof ( tf )  )