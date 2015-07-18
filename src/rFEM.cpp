
#define R_VERSION_

#include "rFEM.hpp"
//#include "IO_handler.hpp"
#include "regressionData.hpp"
#include "mesh_objects.hpp"
#include "mesh.hpp"
#include "finite_element.hpp"
#include "matrix_assembler.hpp"

#include "mixedFERegression.hpp"
#include "evaluator.hpp"


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
	SEXP result;

    if(regressionData.getOrder()==1)
    {
		MeshHandler<1> mesh(Rmesh);
		//std::cout<< "Mesh loaded"<<std::endl;
		MixedFERegression<RegressionData, IntegratorTriangleP2,1> regression(mesh,regressionData);

		regression.smoothLaplace();

		const std::vector<VectorXr>& solution = regression.getSolution();
		const std::vector<Real>& dof = regression.getDOF();
		//Copy result in R memory
		result = PROTECT(allocVector(VECSXP, 2));
		SET_VECTOR_ELT(result, 0, allocMatrix(REALSXP, solution[0].size(), solution.size()));
		SET_VECTOR_ELT(result, 1, allocVector(REALSXP, solution.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(int j = 0; j < solution.size(); j++)
		{
			for(int i = 0; i < solution[0].size(); i++)
				rans[i + solution[0].size()*j] = solution[j][i];
		}

		Real *rans2 = REAL(VECTOR_ELT(result, 1));
		for(int i = 0; i < solution.size(); i++)
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
		result = PROTECT(allocVector(VECSXP, 2));
		SET_VECTOR_ELT(result, 0, allocMatrix(REALSXP, solution[0].size(), solution.size()));
		SET_VECTOR_ELT(result, 1, allocVector(REALSXP, solution.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(int j = 0; j < solution.size(); j++)
		{
			for(int i = 0; i < solution[0].size(); i++)
				rans[i + solution[0].size()*j] = solution[j][i];
		}

		Real *rans2 = REAL(VECTOR_ELT(result, 1));
		for(int i = 0; i < solution.size(); i++)
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
	SEXP result;

    if(regressionData.getOrder()==1)
    {
		MeshHandler<1> mesh(Rmesh);
		//std::cout<< "Mesh loaded"<<std::endl;
		MixedFERegression<RegressionDataElliptic, IntegratorTriangleP2,1> regression(mesh,regressionData);

		regression.smoothEllipticPDE();

		const std::vector<VectorXr>& solution = regression.getSolution();
		const std::vector<Real>& dof = regression.getDOF();
		//Copy result in R memory
		result = PROTECT(allocVector(VECSXP, 2));
		SET_VECTOR_ELT(result, 0, allocMatrix(REALSXP, solution[0].size(), solution.size()));
		SET_VECTOR_ELT(result, 1, allocVector(REALSXP, solution.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(int j = 0; j < solution.size(); j++)
		{
			for(int i = 0; i < solution[0].size(); i++)
				rans[i + solution[0].size()*j] = solution[j][i];
		}

		Real *rans2 = REAL(VECTOR_ELT(result, 1));
		for(int i = 0; i < solution.size(); i++)
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
		result = PROTECT(allocVector(VECSXP, 2));
		SET_VECTOR_ELT(result, 0, allocMatrix(REALSXP, solution[0].size(), solution.size()));
		SET_VECTOR_ELT(result, 1, allocVector(REALSXP, solution.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(int j = 0; j < solution.size(); j++)
		{
			for(int i = 0; i < solution[0].size(); i++)
				rans[i + solution[0].size()*j] = solution[j][i];
		}

		Real *rans2 = REAL(VECTOR_ELT(result, 1));
		for(int i = 0; i < solution.size(); i++)
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
	SEXP result;

    if(regressionData.getOrder()==1)
    {
		MeshHandler<1> mesh(Rmesh);
		//std::cout<< "Mesh loaded"<<std::endl;
		MixedFERegression<RegressionDataEllipticSpaceVarying, IntegratorTriangleP2,1> regression(mesh,regressionData);

		regression.smoothEllipticPDESpaceVarying();

		const std::vector<VectorXr>& solution = regression.getSolution();
		const std::vector<Real>& dof = regression.getDOF();
		//Copy result in R memory
		result = PROTECT(allocVector(VECSXP, 2));
		SET_VECTOR_ELT(result, 0, allocMatrix(REALSXP, solution[0].size(), solution.size()));
		SET_VECTOR_ELT(result, 1, allocVector(REALSXP, solution.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(int j = 0; j < solution.size(); j++)
		{
			for(int i = 0; i < solution[0].size(); i++)
				rans[i + solution[0].size()*j] = solution[j][i];
		}

		Real *rans2 = REAL(VECTOR_ELT(result, 1));
		for(int i = 0; i < solution.size(); i++)
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
		result = PROTECT(allocVector(VECSXP, 2));
		SET_VECTOR_ELT(result, 0, allocMatrix(REALSXP, solution[0].size(), solution.size()));
		SET_VECTOR_ELT(result, 1, allocVector(REALSXP, solution.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(int j = 0; j < solution.size(); j++)
		{
			for(int i = 0; i < solution[0].size(); i++)
				rans[i + solution[0].size()*j] = solution[j][i];
		}

		Real *rans2 = REAL(VECTOR_ELT(result, 1));
		for(int i = 0; i < solution.size(); i++)
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
    	PROTECT(result=allocVector(REALSXP, 2*IntegratorTriangleP2::NNODES*mesh.num_triangles()));

    	FiniteElement<IntegratorTriangleP2,1> fe;
    	for(auto i=0; i<mesh.num_triangles(); i++)
    	{
    		fe.updateElement(mesh.getTriangle(i));
    		for(int l = 0;l < IntegratorTriangleP2::NNODES; l++)
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
    	PROTECT(result=allocVector(REALSXP, 2*IntegratorTriangleP4::NNODES*mesh.num_triangles()));

    	FiniteElement<IntegratorTriangleP4,2> fe;
    	for(auto i=0; i<mesh.num_triangles(); i++)
    	{
    		fe.updateElement(mesh.getTriangle(i));
    		for(int l = 0;l < IntegratorTriangleP4::NNODES; l++)
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

	int n_coef 	= length(Rcoef);
	int n_X 	= length(RX);

    // Cast all computation parameters
    X 			= REAL(RX);
    Y 			= REAL(RY);
    coef 		= REAL(Rcoef);
    order 		= INTEGER(Rorder)[0];
    fast 		= INTEGER(Rfast)[0];

	//std::cout<<"Starting Evaluation"<<std::endl;

    SEXP result;
	PROTECT(result=allocVector(REALSXP, n_X));

    //Set the mesh

    if(order == 1)
    {
		Evaluator<1> evaluator(Rmesh);
		evaluator.eval(X, Y, n_X, coef, order, fast, REAL(result));
	}
	else if(order == 2)
	{
		Evaluator<2> evaluator(Rmesh);
		evaluator.eval(X, Y, n_X, coef, order, fast, REAL(result));
	}


	UNPROTECT(1);
    // result list
    return(result);
}
}



