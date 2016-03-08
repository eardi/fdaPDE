
#define R_VERSION_

#include "fdaPDE.h"
//#include "IO_handler.hpp"
#include "regressionData.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "finite_element.h"
#include "matrix_assembler.h"

#include "mixedFERegression.h"

extern "C" {
//! This function manages the various options for Spatial Regression, Sangalli et al version
/*!
	This function is than called from R code.
	\param Robservations an R-vector containing the values of the observations.
	\param Rdesmat an R-matrix containing the design matrix for the regression.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param Rbindex an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param Rbvalues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
	\return R-vector containg the coefficients of the solution
*/

SEXP regression_Laplace(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rlambda,
				   SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP DOF)
{
    //Set data
	RegressionData regressionData(Rlocations, Robservations, Rorder, Rlambda, Rcovariates, RBCIndices, RBCValues, DOF);

	//std::cout<< "Data loaded"<<std::endl;
	SEXP result = NILSXP;

    if(regressionData.getOrder()==1)
    {
		MeshHandler<1> mesh(Rmesh);
		//std::cout<< "Mesh loaded"<<std::endl;
		MixedFERegression<RegressionData, IntegratorTriangleP2,1> regression(mesh,regressionData);

		regression.smoothLaplace();

		const std::vector<VectorXr>& solution = regression.getSolution();
		const std::vector<Real>& dof = regression.getDOF();
		//Copy result in R memory
		result = PROTECT(Rf_allocVector(VECSXP, 2));
		SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution[0].size(), solution.size()));
		SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, solution.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(UInt j = 0; j < solution.size(); j++)
		{
			for(UInt i = 0; i < solution[0].size(); i++)
				rans[i + solution[0].size()*j] = solution[j][i];
		}

		Real *rans2 = REAL(VECTOR_ELT(result, 1));
		for(UInt i = 0; i < solution.size(); i++)
		{
			rans2[i] = dof[i];
		}
		UNPROTECT(1);

    }
	else if(regressionData.getOrder()==2)
	{
		MeshHandler<2> mesh(Rmesh);
		//std::cout<< "Mesh loaded"<<std::endl;
		MixedFERegression<RegressionData, IntegratorTriangleP4,2> regression(mesh,regressionData);

		regression.smoothLaplace();

		const std::vector<VectorXr>& solution = regression.getSolution();
		const std::vector<Real>& dof = regression.getDOF();
		//Copy result in R memory
		result = PROTECT(Rf_allocVector(VECSXP, 2));
		SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution[0].size(), solution.size()));
		SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, solution.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(UInt j = 0; j < solution.size(); j++)
		{
			for(UInt i = 0; i < solution[0].size(); i++)
				rans[i + solution[0].size()*j] = solution[j][i];
		}

		Real *rans2 = REAL(VECTOR_ELT(result, 1));
		for(UInt i = 0; i < solution.size(); i++)
		{
			rans2[i] = dof[i];
		}
		UNPROTECT(1);
    }

	return(result);
}

SEXP regression_PDE(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc,
				   SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP DOF)
{
    //Set data
	RegressionDataElliptic regressionData(Rlocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Rcovariates, RBCIndices, RBCValues, DOF);

	//std::cout<< "Data loaded"<<std::endl;
	SEXP result = NILSXP;

    if(regressionData.getOrder()==1)
    {
		MeshHandler<1> mesh(Rmesh);
		//std::cout<< "Mesh loaded"<<std::endl;
		MixedFERegression<RegressionDataElliptic, IntegratorTriangleP2,1> regression(mesh,regressionData);

		regression.smoothEllipticPDE();

		const std::vector<VectorXr>& solution = regression.getSolution();
		const std::vector<Real>& dof = regression.getDOF();
		//Copy result in R memory
		result = PROTECT(Rf_allocVector(VECSXP, 2));
		SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution[0].size(), solution.size()));
		SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, solution.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(UInt j = 0; j < solution.size(); j++)
		{
			for(UInt i = 0; i < solution[0].size(); i++)
				rans[i + solution[0].size()*j] = solution[j][i];
		}

		Real *rans2 = REAL(VECTOR_ELT(result, 1));
		for(UInt i = 0; i < solution.size(); i++)
		{
			rans2[i] = dof[i];
		}
		UNPROTECT(1);


    }
	else if(regressionData.getOrder()==2)
	{
		MeshHandler<2> mesh(Rmesh);
		//std::cout<< "Mesh loaded"<<std::endl;
		MixedFERegression<RegressionDataElliptic, IntegratorTriangleP4,2> regression(mesh,regressionData);

		regression.smoothEllipticPDE();

		const std::vector<VectorXr>& solution = regression.getSolution();
		const std::vector<Real>& dof = regression.getDOF();
		//Copy result in R memory
		result = PROTECT(Rf_allocVector(VECSXP, 2));
		SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution[0].size(), solution.size()));
		SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, solution.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(UInt j = 0; j < solution.size(); j++)
		{
			for(UInt i = 0; i < solution[0].size(); i++)
				rans[i + solution[0].size()*j] = solution[j][i];
		}

		Real *rans2 = REAL(VECTOR_ELT(result, 1));
		for(UInt i = 0; i < solution.size(); i++)
		{
			rans2[i] = dof[i];
		}
		UNPROTECT(1);
    }

	return(result);
}



SEXP regression_PDE_space_varying(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru,
				   SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP DOF)
{
    //Set data
	RegressionDataEllipticSpaceVarying regressionData(Rlocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Ru, Rcovariates, RBCIndices, RBCValues, DOF);

	//regressionData.print();

	//std::cout<< "Data Loaded"<<std::endl;
	SEXP result = NILSXP;

    if(regressionData.getOrder()==1)
    {
		MeshHandler<1> mesh(Rmesh);
		//std::cout<< "Mesh loaded"<<std::endl;
		MixedFERegression<RegressionDataEllipticSpaceVarying, IntegratorTriangleP2,1> regression(mesh,regressionData);

		regression.smoothEllipticPDESpaceVarying();

		const std::vector<VectorXr>& solution = regression.getSolution();
		const std::vector<Real>& dof = regression.getDOF();
		//Copy result in R memory
		result = PROTECT(Rf_allocVector(VECSXP, 2));
		SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution[0].size(), solution.size()));
		SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, solution.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(UInt j = 0; j < solution.size(); j++)
		{
			for(UInt i = 0; i < solution[0].size(); i++)
				rans[i + solution[0].size()*j] = solution[j][i];
		}

		Real *rans2 = REAL(VECTOR_ELT(result, 1));
		for(UInt i = 0; i < solution.size(); i++)
		{
			rans2[i] = dof[i];
		}
		UNPROTECT(1);

    }
	else if(regressionData.getOrder()==2)
	{
		MeshHandler<2> mesh(Rmesh);
		//std::cout<< "Mesh loaded"<<std::endl;
		MixedFERegression<RegressionDataEllipticSpaceVarying, IntegratorTriangleP4,2> regression(mesh,regressionData);

		regression.smoothEllipticPDESpaceVarying();

		const std::vector<VectorXr>& solution = regression.getSolution();
		const std::vector<Real>& dof = regression.getDOF();
		//Copy result in R memory
		result = PROTECT(Rf_allocVector(VECSXP, 2));
		SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution[0].size(), solution.size()));
		SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, solution.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(UInt j = 0; j < solution.size(); j++)
		{
			for(UInt i = 0; i < solution[0].size(); i++)
				rans[i + solution[0].size()*j] = solution[j][i];
		}

		Real *rans2 = REAL(VECTOR_ELT(result, 1));
		for(UInt i = 0; i < solution.size(); i++)
		{
			rans2[i] = dof[i];
		}
		UNPROTECT(1);
    }

	return(result);
}

SEXP get_integration_points(SEXP Rmesh, SEXP Rorder)
{
	//Declare pointer to access data from C++

	int order;

	// Cast all computation parameters
    order 		= INTEGER(Rorder)[0];

	//std::cout<<"Computing Locations for Numeric Integration"<<std::endl;

    SEXP result;

    if(order == 1)
    {
    	MeshHandler<1> mesh(Rmesh);
    	PROTECT(result=Rf_allocVector(REALSXP, 2*IntegratorTriangleP2::NNODES*mesh.num_triangles()));

    	FiniteElement<IntegratorTriangleP2,1> fe;
    	for(UInt i=0; i<mesh.num_triangles(); i++)
    	{
    		fe.updateElement(mesh.getTriangle(i));
    		for(UInt l = 0;l < IntegratorTriangleP2::NNODES; l++)
    		{
    			Point p = fe.coorQuadPt(l);
    			REAL(result)[i*IntegratorTriangleP2::NNODES + l] = p[0];
    			REAL(result)[mesh.num_triangles()*IntegratorTriangleP2::NNODES + i*IntegratorTriangleP2::NNODES + l] = p[1];
    		}

    	}
    }
    else if(order == 2)
    {
    	MeshHandler<2> mesh(Rmesh);
    	PROTECT(result=Rf_allocVector(REALSXP, 2*IntegratorTriangleP4::NNODES*mesh.num_triangles()));

    	FiniteElement<IntegratorTriangleP4,2> fe;
    	for(UInt i=0; i<mesh.num_triangles(); i++)
    	{
    		fe.updateElement(mesh.getTriangle(i));
    		for(UInt l = 0;l < IntegratorTriangleP4::NNODES; l++)
    		{
    			Point p = fe.coorQuadPt(l);
    			REAL(result)[i*IntegratorTriangleP4::NNODES + l] = p[0];
    			REAL(result)[mesh.num_triangles()*IntegratorTriangleP4::NNODES + i*IntegratorTriangleP4::NNODES + l] = p[1];
    		}

    	}
    }


	UNPROTECT(1);
    // result list
    return(result);
}


}



