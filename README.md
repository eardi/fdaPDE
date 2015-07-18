Spatial Regression penalized with Partial Differential Equations.
========

This is an R library, implementated in C++, for Spatial Spline Regression Models.
We use an Expression Templates approach for the discretization of the PDE
penalizing the regression.


Further *methodological* details in:

**- [Spatial Spline Regression Models] (http://onlinelibrary.wiley.com/doi/10.1111/rssb.12009/abstract)**

**- [Blood flow velocity field estimation via spatial regression with PDE penalization] (http://mox.polimi.it/it/progetti/pubblicazioni/quaderni/19-2013.pdf)**

Requirements
===================

* A C++11 compiler
* R (>= 1.9.0)
* R Packages: fda, rgl, TriLibrary


Full Installation
=================

1. Clone the git repository with `git clone https://github.com/OldFoxes/SSR.git`
2. Run `R CMD build SSR` where `SSR` is the directory where the repository has been cloned, the command will create the R package `SSR_0.1-2.tar.gz`
3. Open your R session, and run `install.packages("SSR_0.1-2.tar.gz", repos = NULL, type = "source")`
4. For a simple test, run the script in `SSR/R script/SSR.R`. Remember to properly set your working directory!
