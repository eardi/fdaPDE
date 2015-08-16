library(FEMr)
#setwd("~/workspace/RPDE/RScripts")

data(MeuseData)
data(MeuseBorder)

order=2
mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = order)
plot(mesh)

mesh <- refine.MESH.2D(mesh, maximum_area = 5000)
plot(mesh)

data = log(MeuseData[,7])

covariates = matrix(1,nrow=length(data),ncol=2)
covariates[,1] = sqrt(MeuseData[,9])
covariates[,2] = (MeuseData[,8])

#  set up the FEM basis object and plot it

# NO ORDER
basisobj = create.FEM.basis(mesh, order)

#  smooth the data without covariates
lambda = 10^3.5

Rprof("smooth.out", memory.profiling = TRUE)
ouputR<-smooth.FEM.basis(locations = MeuseData[,c(2,3)],
                             observations = data, basisobj = basisobj,
                            lambda = lambda, GCV = FALSE, CPP_CODE = FALSE)
Rprof(NULL)
system.time(smooth.FEM.basis(locations = MeuseData[,c(2,3)],
                             observations = data, basisobj = basisobj,
                             lambda = lambda, GCV = FALSE, CPP_CODE = TRUE))

ZincMeusefd1 = smooth.FEM.basis(locations = MeuseData[,c(2,3)], 
                                observations = data, basisobj = basisobj,
                                  lambda = lambda, GCV = FALSE, CPP_CODE = TRUE)
plot(ZincMeusefd1$felsplobj, num_refinements = 10)

## What if BC Dirichlet

BorderIndices = (1:length(mesh$nodesmarker))[mesh$nodesmarker==1]
Indices= BorderIndices
Values = rep(x = 7, times = length(Indices))
BC = list(Indices = Indices, Values = Values)

output = smooth.FEM.basis(#locations = MeuseData[,c(2,3)], 
                            observations = data, 
                            #covariates = covariates, 
                            basisobj = basisobj,
                            GCV = FALSE,
                            lambda = lambda, BC = BC, CPP_CODE = TRUE)
plot(output$felsplobj, num_refinements = 10)

## Generalized Elliptic PDE (example with preferential smoothing by advection i direction (1,1))

PDE_parameters = list(K = 1*matrix(c(1,0,0,1), nrow = 2), beta = c(0.1,0.1),
                      #c = 50e-4
                      c = 0
                      )
output = smooth.FEM.PDE.basis(observations = data, 
                          basisobj = basisobj, lambda = lambda, PDE_parameters = PDE_parameters, 
                          covariates = NULL,
                          BC = NULL,
                          GCV = FALSE
                          )
plot(output$felsplobj, num_refinements = 10)

