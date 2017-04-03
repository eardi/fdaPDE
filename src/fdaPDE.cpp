
#define R_VERSION_

#include "fdaPDE.h"
//#include "IO_handler.hpp"
#include "regressionData.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "finite_element.h"
#include "matrix_assembler.h"

#include "mixedFERegression.h"


template<typename InputHandler, typename Integrator, UInt ORDER>
SEXP regression_skeleton(InputHandler &regressionData, SEXP Rmesh)
{
	MeshHandler<ORDER> mesh(Rmesh);
	MixedFERegression<InputHandler, Integrator,ORDER> regression(mesh,regressionData);

	regression.apply();

	const std::vector<VectorXr>& solution = regression.getSolution();
	const std::vector<Real>& dof = regression.getDOF();

	//Copy result in R memory
	SEXP result = NILSXP;
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
	return(result);
}

template<typename Integrator, UInt ORDER>
SEXP get_integration_points_skeleton(SEXP Rmesh)
{
	MeshHandler<ORDER> mesh(Rmesh);
	FiniteElement<Integrator,ORDER> fe;

	SEXP result;
	PROTECT(result=Rf_allocVector(REALSXP, 2*Integrator::NNODES*mesh.num_triangles()));
	for(UInt i=0; i<mesh.num_triangles(); i++)
	{
		fe.updateElement(mesh.getTriangle(i));
		for(UInt l = 0;l < Integrator::NNODES; l++)
		{
			Point p = fe.coorQuadPt(l);
			REAL(result)[i*Integrator::NNODES + l] = p[0];
			REAL(result)[mesh.num_triangles()*Integrator::NNODES + i*Integrator::NNODES + l] = p[1];
		}
	}

	UNPROTECT(1);
	return(result);
}

template<typename Integrator, UInt ORDER, typename A>
SEXP get_FEM_Matrix_skeleton(SEXP Rmesh, EOExpr<A> oper)
{
	MeshHandler<ORDER> mesh(Rmesh);

	UInt nnodes=mesh.num_nodes();
	FiniteElement<Integrator, ORDER> fe;

	SpMat AMat;
	Assembler::operKernel(oper, mesh, fe, AMat);
   	//std::cout << AMat;

	//Copy result in R memory
	SEXP result;
	result = PROTECT(Rf_allocVector(VECSXP, 2));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(INTSXP, AMat.nonZeros() , 2));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, AMat.nonZeros()));

	int *rans = INTEGER(VECTOR_ELT(result, 0));
	Real  *rans2 = REAL(VECTOR_ELT(result, 1));
	UInt i = 0;
	for (UInt k=0; k < AMat.outerSize(); ++k)
		{
			for (SpMat::InnerIterator it(AMat,k); it; ++it)
			{
				//std::cout << "(" << it.row() <<","<< it.col() <<","<< it.value() <<")\n";
				rans[i] = 1+it.row();
				rans[i + AMat.nonZeros()] = 1+it.col();
				rans2[i] = it.value();
				i++;
			}
		}
	UNPROTECT(1);
	return(result);
}

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
    //Set input data
	RegressionData regressionData(Rlocations, Robservations, Rorder, Rlambda, Rcovariates, RBCIndices, RBCValues, DOF);

    if(regressionData.getOrder()==1)
    	return(regression_skeleton<RegressionData,IntegratorTriangleP2, 1>(regressionData, Rmesh));
    else if(regressionData.getOrder()==2)
		return(regression_skeleton<RegressionData,IntegratorTriangleP4, 2>(regressionData, Rmesh));
    return(NILSXP);
}

SEXP regression_PDE(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc,
				   SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP DOF)
{
	RegressionDataElliptic regressionData(Rlocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Rcovariates, RBCIndices, RBCValues, DOF);

	if(regressionData.getOrder() == 1)
		return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP2, 1>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 2)
		return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP4, 2>(regressionData, Rmesh));
	return(NILSXP);
}



SEXP regression_PDE_space_varying(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru,
				   SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP DOF)
{
    //Set data
	RegressionDataEllipticSpaceVarying regressionData(Rlocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Ru, Rcovariates, RBCIndices, RBCValues, DOF);

	if(regressionData.getOrder() == 1)
		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP2, 1>(regressionData, Rmesh));
	else if(regressionData.getOrder() == 2)
		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP4, 2>(regressionData, Rmesh));
	return(NILSXP);
}

SEXP get_integration_points(SEXP Rmesh, SEXP Rorder)
{
	//Declare pointer to access data from C++
	int order = INTEGER(Rorder)[0];

    if(order == 1)
    	return(get_integration_points_skeleton<IntegratorTriangleP2, 1>(Rmesh));
    else if(order == 2)
    	return(get_integration_points_skeleton<IntegratorTriangleP4, 2>(Rmesh));
    return(NILSXP);
}

SEXP get_FEM_PDE_Matrix(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc,
				   SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP DOF)
{
	RegressionDataElliptic regressionData(Rlocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Rcovariates, RBCIndices, RBCValues, DOF);

	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
	typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	const Real& c = regressionData.getC();
	const Eigen::Matrix<Real,2,2>& K = regressionData.getK();
	const Eigen::Matrix<Real,2,1>& beta = regressionData.getBeta();

    if(regressionData.getOrder()==1)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	if(regressionData.getOrder()==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	return(NILSXP);
}

SEXP get_FEM_mass_matrix(SEXP Rmesh, SEXP Rorder)
{
	int order = INTEGER(Rorder)[0];
	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);

    if(order==1)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1>(Rmesh, mass));
	if(order==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2>(Rmesh, mass));
	return(NILSXP);
}

SEXP get_FEM_stiff_matrix(SEXP Rmesh, SEXP Rorder)
{
	int order = INTEGER(Rorder)[0];
	typedef EOExpr<Stiff> ETMass;   Stiff EStiff;   ETMass stiff(EStiff);

    if(order==1)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1>(Rmesh, stiff));
	if(order==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2>(Rmesh, stiff));
	return(NILSXP);
}

SEXP get_FEM_PDE_matrix(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc,
				   SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP DOF)
{
	RegressionDataElliptic regressionData(Rlocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Rcovariates, RBCIndices, RBCValues, DOF);

	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
	typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	const Real& c = regressionData.getC();
	const Eigen::Matrix<Real,2,2>& K = regressionData.getK();
	const Eigen::Matrix<Real,2,1>& beta = regressionData.getBeta();

    if(regressionData.getOrder()==1)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	if(regressionData.getOrder()==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	return(NILSXP);
}

SEXP get_FEM_PDE_space_varying_matrix(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru,
		   SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP DOF)
{
	RegressionDataEllipticSpaceVarying regressionData(Rlocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Ru, Rcovariates, RBCIndices, RBCValues, DOF);

	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
	typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	const Reaction& c = regressionData.getC();
	const Diffusivity& K = regressionData.getK();
	const Advection& beta = regressionData.getBeta();

    if(regressionData.getOrder()==1)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	if(regressionData.getOrder()==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	return(NILSXP);
}

}
