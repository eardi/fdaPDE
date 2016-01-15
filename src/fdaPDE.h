#ifndef FDAPDE_H_
#define FDAPDE_H_

// Insert principal libraries
#ifdef R_VERSION_
#define R_NO_REMAP
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h> 
#endif

#include <stdint.h>

#include <cstdlib>
//#include <iomanip>
#include <limits>
#include <vector>

// For debugging purposes
//#include <Eigen/StdVector>
//#include "Eigen/Eigen/Sparse"
//#include "Eigen/Eigen/Dense"
//#define  EIGEN_MPL2_ONLY

//Take the code from the linked RcppEigen
#include <Eigen/StdVector>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#define  EIGEN_MPL2_ONLY

typedef double Real;
typedef int UInt;

typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> MatrixXr;
typedef Eigen::Matrix<Real,Eigen::Dynamic,1> VectorXr;
typedef Eigen::SparseMatrix<Real> SpMat;
typedef Eigen::SparseVector<Real> SpVec;
typedef Eigen::Triplet<Real> coeff;

#endif /* FDAPDE_H_ */
