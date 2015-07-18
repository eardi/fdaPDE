pkgname <- "RPDE"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('RPDE')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("MeuseBorder")
### * MeuseBorder

flush(stderr()); flush(stdout())

### Name: MeuseBorder
### Title: Border of Meuse river data set
### Aliases: MeuseBorder
### Keywords: datasets

### ** Examples

data(MeuseBorder)
summary(MeuseBorder)



cleanEx()
nameEx("MeuseData")
### * MeuseData

flush(stderr()); flush(stdout())

### Name: MeuseData
### Title: Meuse river data set - original, full data set
### Aliases: MeuseData
### Keywords: datasets

### ** Examples

data(MeuseData)
summary(MeuseData)



cleanEx()
nameEx("create.FEM.basis")
### * create.FEM.basis

flush(stderr()); flush(stdout())

### Name: create.FEM.basis
### Title: Creates a Finite Element Method basis
### Aliases: create.FEM.basis

### ** Examples

## Creates an object TRIMESH2D with a concavity and second order nodes
mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                     segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)), order=1)
## Plot it
plot(mesh)                   

# Creates the basis object
basisobj = create.FEM.basis(mesh, order=1)



cleanEx()
nameEx("create.MESH.2D")
### * create.MESH.2D

flush(stderr()); flush(stdout())

### Name: create.MESH.2D
### Title: Creates a TRIMESH2D object
### Aliases: create.MESH.2D

### ** Examples

## Creates an object TRIMESH2D on the convex hull of the specified nodes
mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)))
## Plot it
plot(mesh)

## Creates an object TRIMESH2D with a concavity and second order nodes
mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                     segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)))
## Plot it
plot(mesh)

## Creates an object TRIMESH2D with second order nodes starting from a first order triangulation
## specified by nodes and triangles
mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                     triangles=rbind(c(2,1,3), c(3,5,4), c(5,3,1)), order = 2)
## Plot it
plot(mesh)



cleanEx()
nameEx("refine.MESH.2D")
### * refine.MESH.2D

flush(stderr()); flush(stdout())

### Name: refine.MESH.2D
### Title: Refine the triangulation
### Aliases: refine.MESH.2D

### ** Examples

## Creates an object TRIMESH2D with a concavity and second order nodes
mesh_coarse<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                     segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)))
## Plot it
plot(mesh_coarse)                   
## Refines the the triangulation in specified in the \code{mesh_coarse} object
mesh<-refine.MESH.2D(mesh_coarse,maximum_area = 0.005, delaunay = TRUE)
## Plot the refined mesh
plot(mesh)



cleanEx()
nameEx("smooth.FEM.PDE.SV.basis")
### * smooth.FEM.PDE.SV.basis

flush(stderr()); flush(stdout())

### Name: smooth.FEM.PDE.SV.basis
### Title: Compute a solution for a Spatial Spline problem
### Aliases: smooth.FEM.PDE.SV.basis

### ** Examples

## Creates an object TRIMESH2D with a concavity and second order nodes
mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                     segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)), order=1)
## Plot it
plot(mesh)                   

# Creates the basis object
basisobj = create.FEM.basis(mesh, order=1)



cleanEx()
nameEx("smooth.FEM.PDE.basis")
### * smooth.FEM.PDE.basis

flush(stderr()); flush(stdout())

### Name: smooth.FEM.PDE.basis
### Title: Compute a solution for a Spatial Regression with PDE
###   Penalization
### Aliases: smooth.FEM.PDE.basis

### ** Examples

## Creates an object TRIMESH2D with a concavity and second order nodes
mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                     segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)), order=1)
## Plot it
plot(mesh)                   

# Creates the basis object
basisobj = create.FEM.basis(mesh, order=1)



cleanEx()
nameEx("smooth.FEM.basis")
### * smooth.FEM.basis

flush(stderr()); flush(stdout())

### Name: smooth.FEM.basis
### Title: Compute a solution for a Spatial Spline problem
### Aliases: smooth.FEM.basis

### ** Examples

data(meuseall)
## Creates an object TRIMESH2D with a concavity and second order nodes
mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                     segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)), order=1)
## Plot it
plot(mesh)                   

# Creates the basis object
basisobj = create.FEM.basis(mesh, order=1)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
