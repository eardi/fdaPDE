#ifndef __MIXEDFEREGRESSION_HPP__
#define __MIXEDFEREGRESSION_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "regressionData.h"
#include "solver.h"

//! A LinearSystem class: A class for the linear system construction and resolution.

template<typename InputHandler, typename Integrator, UInt ORDER>
class MixedFERegressionBase{
protected:
	static constexpr Real pruning_coeff = 2.2204e-013;
	static constexpr Real dirichlet_penalization = 10e18;
	const MeshHandler<ORDER> &mesh_;
	const InputHandler& regressionData_;
	std::vector<coeff> tripletsData_;

	//SpMat NWblock_;
	SpMat DMat_;
	SpMat AMat_;
	SpMat MMat_;

	SpMat Psi_;
	MatrixXr P_;

	SpMat coeffmatrix_;       //!A Eigen::VectorXr: Stores the system right hand side.
	VectorXr b_;			  //!A Eigen::VectorXr : Stores the system solution
	std::vector<VectorXr> solution_;
	std::vector<Real> dof_;

	void computeBasisEvaluations();
	void computeProjOnCovMatrix();

	//	|DMat | AMat^T  |
	//	|AMat | MMat	|
	void buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat);
	//! A normal member taking two arguments: Dirichlet Boundary condition
	/*!
	 * This member applies Dirichlet boundary conditions on the linear system with the penalization method.
	  \param bcindex is a const reference to vector<int> : the global indexes of the nodes to which the boundary condition has to be applied.
	  \param bcvalues is a const reference to vector<double> : the values of the boundary conditions relative to bcindex.
	*/
	void applyDirichletBC(const vector<int>& bcindex, const vector<Real>& bcvalues);
	void computeDataMatrix(SpMat& DMat);
	void computeDataMatrixByIndices(SpMat& DMat);
	void computeRightHandData(VectorXr& rightHandData);
	void computeDegreesOfFreedom(UInt output_index);
	 //! A template for the system resolution: SpLu, SpQR, SpCholesky,SpConjGrad
	template<typename P>
	void solve(UInt output_index);

public:
	//!A Constructor.
	MixedFERegressionBase(const MeshHandler<ORDER>& mesh, const InputHandler& regressionData):mesh_(mesh), regressionData_(regressionData){};

	template<typename A>
	void apply(EOExpr<A> oper);

	//! A inline member that returns a VectorXr, returns the whole solution_.
	inline std::vector<VectorXr> const & getSolution() const{return solution_;};
	inline std::vector<Real> const & getDOF() const{return dof_;};
};

template<typename InputHandler, typename Integrator, UInt ORDER>
class MixedFERegression : public MixedFERegressionBase<InputHandler, Integrator, ORDER>
{
public:
	MixedFERegression(const MeshHandler<ORDER>& mesh, const InputHandler& regressionData):MixedFERegressionBase<InputHandler, Integrator, ORDER>(mesh, regressionData){};

	void apply()
	{
		std::cout << "Option not implemented! \n";
	}
};



#include "mixedFERegression_imp.h"

#endif
