#ifndef __IO_HANDLER_HPP__
#define __IO_HANDLER_HPP__

#include "mesh_objects.hpp"
#include "RPDE.hpp"

//!  An IO handler class for objects passed from R
/*!
 * This class, given the data from R, convert them in a C++ format, offering a
 * series of method for their access, so isolating the more the possible the specific
 * code for R/C++ data conversion.
*/
class  IOHandler{
	private:
		
		// Design matrix pointer and dimensions
		VectorXr observations_;
		UInt n_obs_;
		
		//Design matrix
		MatrixXr design_matrix_;
		UInt n_;
		UInt p_;  
		
		//Other parameters
		UInt order_;
		Real lambda_;
		
		Real c_;
		std::vector<Real> beta_;
		Eigen::Matrix <Real,2,2> K_;
		VectorXr u_;
		
		std::vector<Point> locations_;
		std::vector<Real> dirichlet_values_;
		std::vector<UInt> dirichlet_indexes_;
		
		#ifdef R_VERSION_
		void setObservations(SEXP Robservations);
		void setDesignMatrix(SEXP Rdesmat);
		void setLocations(SEXP Rlocations);
		#endif
		
	public:
		//! A basic version of the constructor.
		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Robservations an R-vector containing the values of the observations.
			\param Rdesmat an R-matrix containing the design matrix for the regression.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
			\param Rbindex an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
					the other are automatically considered in Neumann Condition.
			\param Rbvalues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
		*/
		
		IOHandler(VectorXr& observations, MatrixXr& design_matrix, UInt order, Real lambda, std::vector<UInt>& dirichlet_indexes, std::vector<Real>& dirichlet_values):
				observations_(observations), design_matrix_(design_matrix), order_(order), lambda_(lambda)
		{
			dirichlet_indexes_ = dirichlet_indexes;
			dirichlet_values_ = dirichlet_values;
		}
		
		#ifdef R_VERSION_
		IOHandler(SEXP Robservations, SEXP Rdesmat, SEXP Rorder, SEXP Rlambda, SEXP Rbindex, SEXP Rbvalues);
		#endif
		//! A complete version of the constructor.
		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Robservations an R-vector containing the values of the observations.
			\param Rdesmat an R-matrix containing the design matrix for the regression.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
			\param Rbindex an R-integer vector containing the indexes of the nodes the user want to apply a Dirichlet Condition,
					the other are automatically considered in Neumann Condition.
			\param Rbvalues an R-double vector containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
			\param Rc an R-double that contains the coefficient of the REACTION term
			\param Rbeta an R-double 2-dim vector that contains the coefficients for the TRANSPORT coefficients.
			\param RK an R-double 2X2 matrix containing the coefficients for a anisotropic DIFFUSION term.
			\param (UNSUPPORTED put it zero) Ru an R-double vector of length #triangles that contaiins the forcing term integrals.
		*/
		#ifdef R_VERSION_
		IOHandler(SEXP Robservations, SEXP Rlocations, SEXP Rorder, SEXP Rlambda, 
				   SEXP Rbindex, SEXP Rbvalues, SEXP Rc, SEXP Rbeta, SEXP RK, SEXP Ru);
		#endif
		
		void printObservations(std::ostream & out);
		void printDesignMatrix(std::ostream & out);
		void printLocations(std::ostream & out);
		
		//! A method returning a reference to the observations vector
		inline VectorXr const & getObservations() const {return observations_;}
		//! A method returning a reference to the design matrix
		inline MatrixXr const & getDesignMatrix() const {return design_matrix_;}
		//! A method returning the number of observations
		inline UInt const getNumberofObservations() const {return n_obs_;}
		//! A method returning the locations of the observations
		inline std::vector<Point> const & getLocations() const {return locations_;}
		//! A method returning the the penalization term
		inline Real const getLambda() const {return lambda_;}
		//! A method returning the input order
		inline UInt const getOrder() const {return order_;}
		//! A method returning the locations of the observations
		inline Real const getC() const {return c_;}
		//! A method returning the transport coefficient term
		inline std::vector<Real> const & getBeta() const {return beta_;}
		//! A method returning the matrix with unisotropic diffusion coefficients
		inline Eigen::Matrix <Real,2,2> const & getK() const {return K_;}
		//! A method returning the integrals of the forcing term 
		inline VectorXr const & getU() const {return u_;}
		//! A method returning the indexes of the nodes for which is needed to apply Dirichlet Conditions
		inline std::vector<UInt> const & getDirichletIndexes() const {return dirichlet_indexes_;}
		//! A method returning the values to apply for Dirichlet Conditions
		inline std::vector<Real> const & getDirichletValues() const {return dirichlet_values_;}
};

#include "IO_handler_imp_old.hpp"

#endif
