/*
 * FEMeval.cpp
 *
 *  Created on: Aug 16, 2015
 *      Author: eardi
 */


#define R_VERSION_

#include "fdaPDE.h"
//#include "IO_handler.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "evaluator.h"

extern "C" {
//! This function manages the various option for the solution evaluation.
/*!
	This function is than one called from R code.
	Call's the walking algoritm for efficient point location inside the mesh.

	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param RX an R-vector containing the x coordinates of the points to be evaluated
	\param RY an R-vector containing the y coordinates of the points to be evaluated
	\param Rcoef an R-vector the coeficients of the solution
	\param Rorder an R integer containg the order of the solution
	\param Rfast an R integer 0 for Naive location algorithm, 1 for Walking Algorithm (can miss location for non convex meshes)
*/

SEXP eval_FEM_fd(SEXP Rmesh, SEXP RX, SEXP RY, SEXP Rcoef, SEXP Rorder, SEXP Rfast)
{
	//Declare pointer to access data from C++

    double *X, *Y, *coef;
	int order;
	bool fast;

	//int n_coef 	= Rf_length(Rcoef);
	int n_X 	= Rf_length(RX);

    // Cast all computation parameters
    X 			= REAL(RX);
    Y 			= REAL(RY);
    coef 		= REAL(Rcoef);
    order 		= INTEGER(Rorder)[0];
    fast 		= INTEGER(Rfast)[0];

    SEXP result;
	PROTECT(result=Rf_allocVector(REALSXP, n_X));
	std::vector<bool> isinside(n_X);
    //Set the mesh
	//std::cout<<"Length "<<n_X<<"--X0 "<<X[0]<<"--Y0 "<<Y[0];
    if(order == 1)
    {
    	MeshHandler<1> mesh(Rmesh);
		Evaluator<1> evaluator(mesh);
		//std::cout<<"Starting evaluation from FEMeval \n";
		evaluator.eval(X, Y, n_X, coef, order, fast, REAL(result), isinside);
	}
	else if(order == 2)
	{
    	MeshHandler<2> mesh(Rmesh);
    	Evaluator<2> evaluator(mesh);
		evaluator.eval(X, Y, n_X, coef, order, fast, REAL(result), isinside);
	}

    for (int i=0; i<n_X;++i)
    {
    	if(!(isinside[i]))
    	{
    		REAL(result)[i]=NA_REAL;
    	}

    }

	UNPROTECT(1);
    // result list
    return(result);
}
}




