#ifndef __EVALUATOR_HPP__
#define __EVALUATOR_HPP__

#include <iostream>

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"

//!  A class for the evaluation of the solution, given the coefficients of the global bases
/*!
 * This class, given a vector of coordinates evaluates the solution, fist locating the point 
 * and secondly evaluating it respect to the coefficients of the bases
 * It dependes on a template parameter that speciefies the Order of the initializing mesh
*/
template <UInt ORDER>
class Evaluator
{
	public:
		//! A constructor. It initializes the constructor given a mesh object.
		Evaluator(const MeshHandler<ORDER>& mesh): 
								mesh_(mesh){};
		
		//! A member that computes the evaluation of a Point in a mesh, given the bases' coefficients.
		/*!
		\param X a pointer to the x coordinates to evaluate.
		\param Y a pointer to the y coordinates to evaluate.
		\param length a unsigned integer containing the number of points to evaluate.
		\param coef a pointer to the vector of coefficients of the solution, the value in position i 
		is associated to the basis \phi(i)
		\param order a unsigned integer that specifies the order of the solution (1 or 2)
		\param fast a boolean that specifies if the algorithm is completely based on the walking
				algorithm (can miss locations in case of non convex structures)
		\param result a double pointer to an already allocated memory space, where the evaluations
		will be stored
		*/
		void eval(Real* X, Real *Y, UInt length, const Real *coef, UInt order, bool redundancy, Real* result, std::vector<bool>& isinside);
		
	private:
		const MeshHandler<ORDER> &mesh_;

};

#include "evaluator_imp.h"


#endif
