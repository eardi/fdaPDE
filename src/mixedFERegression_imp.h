#ifndef __MIXEDFEREGRESSION_IMP_HPP__
#define __MIXEDFEREGRESSION_IMP_HPP__

#include <iostream>

template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegressionBase<InputHandler,Integrator,ORDER>::computeBasisEvaluations(){

	//std::cout<<"Data Matrix Computation by Basis Evaluation.."<<std::endl;
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();
	Real eps = 2.2204e-016,
		 tolerance = 100 * eps;

	Psi_.resize(nlocations, nnodes);
	//Psi_.reserve(Eigen::VectorXi::Constant(nlocations,ORDER*3));

	Triangle<ORDER*3> tri_activated;
	Eigen::Matrix<Real,ORDER * 3,1> coefficients;

	Real evaluator;
	for(UInt i=0; i<nlocations;i++)
	{
		tri_activated = mesh_.findLocationNaive(regressionData_.getLocations()[i]);
		if(tri_activated.getId() == Identifier::NVAL)
		{
			#ifdef R_VERSION_
			Rprintf("WARNING: Observation %d is not in the domain, remove point and re-perform smoothing\n", i+1);
			#else
			std::cout << "WARNING: Observation " << i+1 <<" is not in the domain\n";
			#endif
		}else
		{
			for(UInt node = 0; node < ORDER*3 ; ++node)
			{
				coefficients = Eigen::Matrix<Real,ORDER * 3,1>::Zero();
				coefficients(node) = 1; //Activates only current base
				evaluator = evaluate_point<ORDER>(tri_activated, regressionData_.getLocations()[i], coefficients);
				Psi_.insert(i, tri_activated[node].getId()) = evaluator;
			}
		}
	}

	Psi_.prune(tolerance);
	Psi_.makeCompressed();
}

template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegressionBase<InputHandler,Integrator,ORDER>::computeProjOnCovMatrix()
{
	//std::cout<<"Computing Projection Matrix"<<std::endl;
	UInt nlocations = regressionData_.getNumberofObservations();

	//Set covariate matrix
	MatrixXr W(this->regressionData_.getCovariates());
	if(regressionData_.isLocationsByNodes())
	{
		MatrixXr W_reduced(regressionData_.getNumberofObservations(), W.cols());
		for (auto i=0; i<nlocations;++i)
		{
			auto index_i = regressionData_.getObservationsIndices()[i];
			for (auto j=0; j<W.cols();++j)
			{
				W_reduced(i,j) = W(index_i,j);
			}
		}
		W = W_reduced;
	}

	//std::cout << W << std::endl;
	MatrixXr WTW(W.transpose()*W);
	P_.resize(W.rows(),W.rows());
	P_ = -W*WTW.ldlt().solve(W.transpose());
	for (int i=0; i<P_.rows();++i)
	{
		P_(i,i) += 1;
	}
	//std::cout << P_ << std::endl;
}


template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegressionBase<InputHandler,Integrator,ORDER>::buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat)
{
	UInt nnodes = mesh_.num_nodes();

	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat.nonZeros() + 2*AMat.nonZeros() + MMat.nonZeros());

	for (int k=0; k<DMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(DMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
	  }
	for (int k=0; k<MMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(MMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes,it.value()));
	  }
	for (int k=0; k<AMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(AMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.col(), it.row()+nnodes,it.value()));
	  }
	for (int k=0; k<AMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(AMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
	  }

	coeffmatrix_.setZero();
	coeffmatrix_.resize(2*nnodes,2*nnodes);
	coeffmatrix_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	coeffmatrix_.makeCompressed();
}

template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegressionBase<InputHandler,Integrator,ORDER>::applyDirichletBC(const vector<int>& bcindex, const vector<Real>& bcvalues)
{
	UInt nnodes = mesh_.num_nodes();

	const std::vector<UInt>& bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real>& bc_values = regressionData_.getDirichletValues();
	UInt nbc_indices = bc_indices.size();

	UInt id1,id3;
	for( auto i=0; i<nbc_indices; i++)
	 {
		id1=bcindex[i];
		id3=id1+nnodes;

		coeffmatrix_.coeffRef(id1,id1) = dirichlet_penalization;
		coeffmatrix_.coeffRef(id3,id3) = dirichlet_penalization;


		b_(id1) += bc_values[i]*dirichlet_penalization;
		b_(id3) = 0;
	 }
	coeffmatrix_.makeCompressed();
}

//construct NW block of the system matrix when basis evaluation is necessary
//!! Depends on computeBasisEvaluations and computeProjOnCovMatrix
template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegressionBase<InputHandler,Integrator,ORDER>::computeDataMatrix(SpMat& DMat)
{
		UInt nnodes = mesh_.num_nodes();
		DMat.resize(nnodes,nnodes);

		if (regressionData_.getCovariates().rows() == 0)
			DMat = Psi_.transpose()*Psi_;
		else
			DMat = (SpMat(Psi_.transpose())*P_*Psi_).sparseView();
}

//construct NW block of the system matrix in Ramsay when locations of observations are
//a subset of the meshe's nodes
//!! Depends on computeBasisEvaluations and computeProjOnCovMatrix
template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegressionBase<InputHandler,Integrator,ORDER>::computeDataMatrixByIndices(SpMat& DMat)
{
		UInt nnodes = mesh_.num_nodes();
		UInt nlocations = regressionData_.getNumberofObservations();

		DMat.resize(nnodes,nnodes);

		if (regressionData_.getCovariates().rows() == 0)
		{
			DMat.reserve(1);
			for (auto i = 0; i<nlocations; ++i)
			{
				auto index = regressionData_.getObservationsIndices()[i];
				DMat.insert(index,index) = 1;
			}
		}
		else
		{
			//May be inefficient
			for (auto i = 0; i<nlocations; ++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				for (auto j = 0; j<nlocations; ++j)
				{
					auto index_j = regressionData_.getObservationsIndices()[j];
					DMat.insert(index_i,index_j) = P_(i,j);
				}
			}
		}
}


template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegressionBase<InputHandler,Integrator,ORDER>::computeRightHandData(VectorXr& rightHandData)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();
	//rightHandData.resize(nnodes);
	rightHandData = VectorXr::Zero(nnodes);

	if(regressionData_.getCovariates().rows() == 0)
	{
		if(regressionData_.isLocationsByNodes())
		{

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				rightHandData(index_i) = regressionData_.getObservations()[i];
			}
		}
		else
		{
			rightHandData=Psi_.transpose()*regressionData_.getObservations();
		}
	}
	else
	{
		if(regressionData_.isLocationsByNodes())
		{
			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				rightHandData(index_i) = (P_).row(i) * regressionData_.getObservations();
			}
		}
		else
		{
			rightHandData=Psi_.transpose()*P_*regressionData_.getObservations();
		}
	}
}


template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegressionBase<InputHandler,Integrator,ORDER>::computeDegreesOfFreedom(UInt output_index)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();

	Eigen::SparseLU<SpMat> solver;
	solver.compute(coeffmatrix_);
	SpMat I(coeffmatrix_.rows(),coeffmatrix_.cols());
	I.setIdentity();
	SpMat coeff_inv = solver.solve(I);


	Real degrees=0;

	if(regressionData_.getCovariates().rows() == 0)
	{
		if(regressionData_.isLocationsByNodes())
		{
			VectorXr d = coeff_inv.diagonal();

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				degrees+=d(index_i);
			}
		}
		else
		{
			MatrixXr An(coeff_inv.topLeftCorner(nnodes, nnodes));
			MatrixXr S = Psi_*An*Psi_.transpose();
			for (auto i=0; i<nlocations;++i)
			{
				degrees+=S(i,i);
			}
		}
	}
	else
	{	//Important
		degrees = regressionData_.getCovariates().cols();

		if(regressionData_.isLocationsByNodes())
		{
			MatrixXr S(coeff_inv);

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				for (auto j=0; j<nlocations;++j)
				{
					auto index_j = regressionData_.getObservationsIndices()[j];

					degrees+=S(index_i,index_j)*P_(j,i);
				}
			}

		}
		else
		{
			MatrixXr An(coeff_inv.topLeftCorner(nnodes, nnodes));
			MatrixXr S = Psi_*An*Psi_.transpose()*P_;
			for (auto i=0; i<nlocations;++i)
			{
				degrees+=S(i,i);
			}
		}
	}


	//std::cout<<"TRACE "<<degrees<<std::endl;

	dof_[output_index] = degrees;
}

template<typename InputHandler, typename Integrator, UInt ORDER>
template<typename A>
void MixedFERegressionBase<InputHandler,Integrator,ORDER>::apply(EOExpr<A> oper)
{
	UInt nnodes=mesh_.num_nodes();
	FiniteElement<Integrator, ORDER> fe;

	//std::cout << "HERE1" << std::endl;
	if(!regressionData_.isLocationsByNodes())
		computeBasisEvaluations();
	if(!(regressionData_.getCovariates().rows() == 0))
		computeProjOnCovMatrix();
	if(!regressionData_.isLocationsByNodes())
		computeDataMatrix(DMat_);
	else
		computeDataMatrixByIndices(DMat_);

	//std::cout << "HERE2" << std::endl;
	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	Assembler::operKernel(oper, mesh_, fe, AMat_);
	Assembler::operKernel(mass, mesh_, fe, MMat_);
	//std::cout << "HERE3" << std::endl;
	VectorXr rightHandData;
	computeRightHandData(rightHandData);
	b_ = VectorXr::Zero(2*nnodes);
	b_.topRows(nnodes)=rightHandData;
	//std::cout << "HERE4" << std::endl;
	solution_.resize(regressionData_.getLambda().size());
	dof_.resize(regressionData_.getLambda().size());

	for(UInt i = 0; i<regressionData_.getLambda().size(); ++i)
	{
		Real lambda = regressionData_.getLambda()[i];
		SpMat AMat_lambda = (-lambda)*AMat_;
		SpMat MMat_lambda = (-lambda)*MMat_;
		this->buildCoeffMatrix(DMat_, AMat_lambda, MMat_lambda);
		//std::cout << coeffmatrix_ << std::endl;
		//Applying boundary conditions if necessary
		if(regressionData_.getDirichletIndices().size() != 0)
			applyDirichletBC(regressionData_.getDirichletIndices(), regressionData_.getDirichletValues());

		this-> template solve<SpLU>(i);
		if(regressionData_.computeDOF())
			computeDegreesOfFreedom(i);
		else
			dof_[i] = -1;
	}
	//std::cout << "HERE5" << std::endl;
}

//solve sparse system with P method

template<typename InputHandler, typename Integrator, UInt ORDER>
template <typename P>
void MixedFERegressionBase<InputHandler,Integrator,ORDER>::solve(UInt output_index)
{
	//std::cout<<this->coeffmatrix_;
	this->solution_[output_index].resize(this->coeffmatrix_.rows());
	P::solve(this->coeffmatrix_,this->b_,this->solution_[output_index]);
}

template<typename Integrator, UInt ORDER>
class MixedFERegression<RegressionData, Integrator, ORDER> : public MixedFERegressionBase<RegressionData, Integrator, ORDER>
{
public:
	MixedFERegression(const MeshHandler<ORDER>& mesh, const RegressionData& regressionData):MixedFERegressionBase<RegressionData, Integrator, ORDER>(mesh, regressionData){};

	void apply()
	{
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	    MixedFERegressionBase<RegressionData, Integrator, ORDER>::apply(stiff);
	}
};

template<typename Integrator, UInt ORDER>
class MixedFERegression<RegressionDataElliptic, Integrator, ORDER> : public MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER>
{
public:
	MixedFERegression(const MeshHandler<ORDER>& mesh, const RegressionDataElliptic& regressionData):MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER>(mesh, regressionData){};

	void apply()
	{
		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	    const Real& c = MixedFERegression::regressionData_.getC();
	    const Eigen::Matrix<Real,2,2>& K = MixedFERegression::regressionData_.getK();
	    const Eigen::Matrix<Real,2,1>& beta = MixedFERegression::regressionData_.getBeta();

	    MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER>::apply(c*mass+stiff[K]+dot(beta,grad));
	}
};

template<typename Integrator, UInt ORDER>
class MixedFERegression<RegressionDataEllipticSpaceVarying, Integrator, ORDER> : public MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER>
{
public:
	MixedFERegression(const MeshHandler<ORDER>& mesh, const RegressionDataEllipticSpaceVarying& regressionData):MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER>(mesh, regressionData){};

	void apply()
	{
		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

		const Reaction& c = MixedFERegression::regressionData_.getC();
		const Diffusivity& K = MixedFERegression::regressionData_.getK();
		const Advection& beta = MixedFERegression::regressionData_.getBeta();

		MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER>::apply(c*mass+stiff[K]+dot(beta,grad));
	}
};

#endif
