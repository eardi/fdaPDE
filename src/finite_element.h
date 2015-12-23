#ifndef __FINITE_ELEMENT_HPP__
#define __FINITE_ELEMENT_HPP__

#include "fdaPDE.h"
#include "integration.h"
#include "mesh_objects.h"

//!  This class implements all properties of a Triangular Finite Element
/*!
 * This class is the most important one of the entire code
 * and implements everything needed by a triangular finite elemnt
 * 
 * It takes as a template parameter a class that implements the mathod used
 * for determining the mass, stiff and grad matrices
*/
template <class Integrator ,UInt ORDER>
class FiniteElement{
private:
	Triangle<ORDER*3> reference_;
	Triangle<ORDER*3> t_;
	Eigen::Matrix<Real,3*ORDER, Integrator::NNODES> phiMapMaster_;
	//Numero basi locali x Num coordinate x numero nodi integrazione
	Eigen::Matrix<Real,3*ORDER, Integrator::NNODES*2> phiDerMapMaster_;
	Eigen::Matrix<Real,3*ORDER, Integrator::NNODES*2> invTrJPhiDerMapMaster_;
	
	void setPhiMaster();
	void setPhiDerMaster();
	void setInvTrJPhiDerMaster();

public:

	//! This is an empty constructor
    /*!
        For efficiency and Expression Templates organization of the
        code, the use of this class is based on the updateElement class
    */
	FiniteElement();
	
	//! A member updating the Finite Element properties
    /*!
      \param t a triangle from which to update the finite element properties
    */
	void updateElement(Triangle<ORDER*3> t);
	
	Real getAreaReference()
	{
		return reference_.getArea();
	}

	Real getDet()
	{
		return t_.getDetJ();
	}

	Point coorQuadPt(UInt iq)
	{
		return Point(t_.getM_J()(0,0)*Integrator::NODES[iq][0] + t_.getM_J()(0,1)*Integrator::NODES[iq][1] + t_[0][0],
				t_.getM_J()(1,0)*Integrator::NODES[iq][0] + t_.getM_J()(1,1)*Integrator::NODES[iq][1] + t_[0][1]);
	}
	
	UInt getGlobalIndex(UInt iq)
	{
		return Integrator::NNODES * t_.getId() + iq;
	}

	//Returns \hat{phi}
	Real phiMaster(UInt i, UInt iq) const;

	//Returns \nabla \hat{phi}
	Real phiDerMaster(UInt i, UInt ic, UInt iq) const;

	//Returns J^{-1} \nabla \hat{phi}
	Real invTrJPhiDerMaster(UInt i, UInt ic, UInt iq) const;
};

#include "finite_element_imp.h"

#endif
