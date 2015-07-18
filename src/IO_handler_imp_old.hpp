#ifndef __IO_HANDLER_IMP_HPP__
#define __IO_HANDLER_IMP_HPP__


#ifdef R_VERSION_
IOHandler::
IOHandler(SEXP Robservations, SEXP Rdesmat, SEXP Rorder, SEXP Rlambda, SEXP Rbindex, SEXP Rbvalues)
{
	setDesignMatrix(Rdesmat);
	setObservations(Robservations);
	
	order_ =  INTEGER(Rorder)[0];
    lambda_ =  REAL(Rlambda)[0];
	
	UInt length_indexes = length(Rbindex);
    for (UInt i = 0; i<length_indexes; ++i)  dirichlet_indexes_.push_back(INTEGER(Rbindex)[i]);
    for (UInt i = 0; i<length_indexes; ++i)  dirichlet_values_.push_back(REAL(Rbvalues)[i]);
	
}

IOHandler::
IOHandler(SEXP Robservations, SEXP Rlocations, SEXP Rorder, SEXP Rlambda, 
				   SEXP Rbindex, SEXP Rbvalues, SEXP Rc, SEXP Rbeta, SEXP RK, SEXP Ru)
{
	setLocations(Rlocations);
	setObservations(Robservations);
	
	order_ =  INTEGER(Rorder)[0];
    lambda_ =  REAL(Rlambda)[0];
    c_ = REAL(Rc)[0];
    
    UInt length_indexes = length(Rbindex);
    for (UInt i = 0; i<length_indexes; ++i)  dirichlet_indexes_.push_back(INTEGER(Rbindex)[i]);
    for (UInt i = 0; i<length_indexes; ++i)  dirichlet_values_.push_back(REAL(Rbvalues)[i]);
    
    beta_.push_back(REAL(Rbeta)[0]);
    beta_.push_back(REAL(Rbeta)[1]);
    
    K_ << REAL(RK)[0] , REAL(RK)[2], REAL(RK)[1], REAL(RK)[3];
    
    UInt length_u = length(Ru);
    u_.resize(length_u);
    
    for (UInt i = 0; i<length_u; ++i)		u_[i] = REAL(Ru)[i];
	
}

void IOHandler::setObservations(SEXP Robservations)
{	
	n_obs_ = INTEGER(getAttrib(Robservations, R_DimSymbol))[0];
	observations_.resize(n_obs_);
	
	
	for(auto i=0;i<n_obs_;++i)
	{
		observations_(i) = REAL(Robservations)[i];
	}
	
}

void IOHandler::setDesignMatrix(SEXP Rdesmat)
{
	n_ = INTEGER(getAttrib(Rdesmat, R_DimSymbol))[0];
	p_ = INTEGER(getAttrib(Rdesmat, R_DimSymbol))[1];
	
	design_matrix_.resize(n_, p_);		
	
	for(auto i=0; i<n_; ++i)
	{
		for(auto j=0; j<p_ ; ++j)
		{
			design_matrix_(i,j)=REAL(Rdesmat)[i+ n_*j];
		}
	}
}

void IOHandler::setLocations(SEXP Rlocations)
{
	n_ = INTEGER(getAttrib(Rlocations, R_DimSymbol))[0];
	
	for(auto i=0; i<n_; ++i)
	{
		locations_.emplace_back(REAL(Rlocations)[i+ n_*0],REAL(Rlocations)[i+ n_*1]);
	}
}
#endif

void IOHandler::printObservations(std::ostream & out)
{
	
	for(auto i=0;i<observations_.size(); i++)
	{
		out<<i<<"\t"<<observations_(i)<<std::endl;
	}
}

void IOHandler::printDesignMatrix(std::ostream & out)
{
	
	for(auto i=0;i<design_matrix_.rows(); i++)
	{
		for(auto j=0; j<design_matrix_.cols(); j++)
		{
			out<<design_matrix_(i,j)<<"\t";
		}
		out<<std::endl;
	}
}

void IOHandler::printLocations(std::ostream & out)
{
	
	for(std::vector<Point>::size_type i=0;i<locations_.size(); i++)
	{
		locations_[i].print(out);
		std::cout<<std::endl;
	}
}

#endif
