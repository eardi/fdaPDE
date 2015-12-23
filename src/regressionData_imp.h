#ifndef __REGRESSIONDATA_IMP_HPP__
#define __REGRESSIONDATA_IMP_HPP__

RegressionData::RegressionData(std::vector<Point>& locations, VectorXr& observations, UInt order, std::vector<Real> lambda, MatrixXr& covariates , std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF):
					locations_(locations), observations_(observations), covariates_(covariates), order_(order), lambda_(lambda),
					bc_values_(bc_values), bc_indices_(bc_indices), DOF_(DOF)
{
	if(locations_.size()==0)
	{
		locations_by_nodes_= true;
		for(int i = 0; i<observations_.size();++i) observations_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_= false;
	}
}

RegressionDataElliptic::RegressionDataElliptic(std::vector<Point>& locations, VectorXr& observations, UInt order, std::vector<Real> lambda, Eigen::Matrix<Real,2,2>& K,	Eigen::Matrix<Real,2,1>& beta,
		Real c, MatrixXr& covariates , std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF):
		 RegressionData(locations, observations, order, lambda, covariates , bc_indices, bc_values, DOF), K_(K), beta_(beta), c_(c)
{;}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(std::vector<Point>& locations, VectorXr& observations, UInt order, std::vector<Real> lambda, const std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > >& K,	const std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > >& beta,
		const std::vector<Real>& c, const std::vector<Real>& u, MatrixXr& covariates , std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF):
		RegressionData(locations, observations, order, lambda, covariates , bc_indices, bc_values, DOF), K_(K), beta_(beta), c_(c), u_(u)
{;}


#ifdef R_VERSION_
RegressionData::RegressionData(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP Rcovariates,
			   SEXP RBCIndices, SEXP RBCValues, SEXP DOF)
{
	setLocations(Rlocations);
	//std::cout<< "Locations set"<<std::endl;
	setObservations(Robservations);
	//std::cout<< "Observations set"<<std::endl;
	setCovariates(Rcovariates);
	//std::cout<< "Covariates set"<<std::endl;

	order_ =  INTEGER(Rorder)[0];
	DOF_ = INTEGER(DOF)[0];

	UInt length_indexes = Rf_length(RBCIndices);
    //for (UInt i = 0; i<length_indexes; ++i)  bc_indices_.push_back(INTEGER(RBCIndices)[i]);
    //for (UInt i = 0; i<length_indexes; ++i)  bc_values_.push_back(REAL(RBCValues)[i]);
	bc_indices_.assign(INTEGER(RBCIndices), INTEGER(RBCIndices) +  length_indexes);
	bc_values_.assign(REAL(RBCValues),REAL(RBCValues) + Rf_length(RBCIndices));

    UInt length_lambda = Rf_length(Rlambda);
    for (UInt i = 0; i<length_lambda; ++i)  lambda_.push_back(REAL(Rlambda)[i]);

}

RegressionDataElliptic::RegressionDataElliptic(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta,
				 SEXP Rc, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP DOF):
	RegressionData(Rlocations, Robservations, Rorder, Rlambda, Rcovariates,
					 			   RBCIndices, RBCValues, DOF)
{
	K_.resize(2, 2);
	for(auto i=0; i<2; ++i)
	{
		for(auto j=0; j<2 ; ++j)
		{
			K_(i,j)=REAL(RK)[i+ 2*j];
		}
	}

	beta_.resize(2);
	for(auto i=0; i<2 ; ++i)
	{
		beta_(i)=REAL(Rbeta)[i];
	}

	c_ =  REAL(Rc)[0];
}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta,
				 SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP DOF):
					 RegressionData(Rlocations, Robservations, Rorder, Rlambda, Rcovariates, RBCIndices, RBCValues, DOF),
					 K_(RK), beta_(Rbeta), c_(Rc), u_(Ru)
{;}

void RegressionDataEllipticSpaceVarying::print(std::ostream & out) const
{
	for (auto i=0;i<18;i++)
		out<<K_(i);
	for (auto i=0;i<18;i++)
		out<<beta_(i);
	for (auto i=0;i<18;i++)
		out<<c_(i);
}

void RegressionData::setObservations(SEXP Robservations)
{
	UInt n_obs_ = Rf_length(Robservations);
	observations_.resize(n_obs_);
	observations_indices_.reserve(n_obs_);

	UInt count = 0;
	if(locations_.size() == 0)
	{
		locations_by_nodes_ = true;
		for(auto i=0;i<n_obs_;++i)
		{
			if(!ISNA(REAL(Robservations)[i]))
			{
				observations_[count] = REAL(Robservations)[i];
				count++;
				observations_indices_.push_back(i);
			}
		}
		observations_.conservativeResize(count, Eigen::NoChange);
	}
	else
	{
		locations_by_nodes_ = false;
		for(auto i=0;i<n_obs_;++i)
		{
			observations_[i] = REAL(Robservations)[i];
		}
	}

	//std::cout<<"Observations #"<<observations_.size()<<std::endl<<observations_<<std::endl;
	//for(auto i=0;i<observations_indices_.size();++i)	std::cout<<observations_indices_[i]<<std::endl;
}

void RegressionData::setCovariates(SEXP Rcovariates)
{


	n_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[0];
	p_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[1];

	covariates_.resize(n_, p_);

	for(auto i=0; i<n_; ++i)
	{
		for(auto j=0; j<p_ ; ++j)
		{
			covariates_(i,j)=REAL(Rcovariates)[i+ n_*j];
		}
	}
}

void RegressionData::setLocations(SEXP Rlocations)
{
	n_ = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];

	for(auto i=0; i<n_; ++i)
	{
		locations_.emplace_back(REAL(Rlocations)[i+ n_*0],REAL(Rlocations)[i+ n_*1]);
	}
}

#endif

void RegressionData::printObservations(std::ostream & out) const
{

	for(auto i=0;i<observations_.size(); i++)
	{
		out<<i<<"\t"<<observations_(i)<<std::endl;
	}
}

void RegressionData::printCovariates(std::ostream & out) const
{

	for(auto i=0;i<covariates_.rows(); i++)
	{
		for(auto j=0; j<covariates_.cols(); j++)
		{
			out<<covariates_(i,j)<<"\t";
		}
		out<<std::endl;
	}
}

void RegressionData::printLocations(std::ostream & out) const
{

	for(std::vector<Point>::size_type i=0;i<locations_.size(); i++)
	{
		locations_[i].print(out);
		//std::cout<<std::endl;
	}
}

#endif
