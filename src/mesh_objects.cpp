/*
 * mesh_object.cpp
 *
 *  Created on: Aug 16, 2015
 *      Author: eardi
 */

#include "mesh_objects.h"


const UInt Identifier::NVAL=std::numeric_limits<UInt>::max();


void Point::print(std::ostream & out) const
{
	out<<"Point -"<< id_ <<"- "<<"("<<coord_[0]<<","<<coord_[1]<<")"<<std::endl<<"------"<<std::endl;
}

void Edge::print(std::ostream & out) const
{
	out<<"Edge -"<< id_ <<"- "<<"("<<points_[0].getId()<<","<<points_[1].getId()<<")"<<std::endl;
}

//template <UInt ORDER>
//Real evaluate_point(const Triangle<3*ORDER>& t, const Point& point, const Eigen::Matrix<Real,3*ORDER,1>& coefficients)
//{
//	//std::cerr<< "TRYING TO EVALUATE ORDER NOT IMPLEMENTED" << std::endl;
//	return 0;
//}
//
//template <>
//Real evaluate_point<1>(const Triangle<3>& t, const Point& point, const Eigen::Matrix<Real,3,1>& coefficients)
//{
//	Eigen::Matrix<Real,3,1> bary_coeff = t.getBaryCoordinates(point);
//	//std::cout<< "B-coord: "<<bary_coeff<<std::endl;
//
//	return(coefficients.dot(bary_coeff));
//}
//
//template <>
//Real evaluate_point<2>(const Triangle<6>& t, const Point& point, const Eigen::Matrix<Real,6,1>& coefficients)
//{
//	Eigen::Matrix<Real,3,1> bary_coeff = t.getBaryCoordinates(point);
//	return( coefficients[0]*(2*bary_coeff[0]*bary_coeff[0]- bary_coeff[0]) +
//            coefficients[1]*(2*bary_coeff[1]*bary_coeff[1] - bary_coeff[1]) +
//            coefficients[2]*(2*bary_coeff[2]*bary_coeff[2] - bary_coeff[2]) +
//            coefficients[3]*(4*bary_coeff[1]* bary_coeff[2])    +
//            coefficients[4]*(4*bary_coeff[2]* bary_coeff[0])    +
//            coefficients[5]*(4*bary_coeff[0]* bary_coeff[1]) );
//
//}
//
//template <UInt ORDER>
//Eigen::Matrix<Real,2,1> evaluate_der_point(const Triangle<3*ORDER>& t, const Point& point, const Eigen::Matrix<Real,3*ORDER,1>& coefficients)
//{
//	//std::cerr<< "TRYING TO EVALUATE ORDER NOT IMPLEMENTED" << std::endl;
//	Eigen::Matrix<Real,2,1> null;
//	return(null);
//}
//
//template <>
//Eigen::Matrix<Real,2,1> evaluate_der_point<1>(const Triangle<3>& t, const Point& point, const Eigen::Matrix<Real,3,1>& coefficients)
//{
//	Eigen::Matrix<Real,2,3> B1;
//	B1 << t[1][1] - t[2][1], t[2][1] - t[0][1], t[0][1] - t[1][1],
//		t[2][0] - t[1][0], t[0][0] - t[2][0], t[1][0] - t[0][0];
//	B1 = B1 / (2 * t.getArea());
//
//	return(B1*coefficients);
//
//}
//
//template <>
//Eigen::Matrix<Real,2,1> evaluate_der_point<2>(const Triangle<6>& t, const Point& point, const Eigen::Matrix<Real,6,1>& coefficients)
//{
//	Eigen::Matrix<Real,3,1> L = t.getBaryCoordinates(point);
//	Eigen::Matrix<Real,2,3> B1;
//	B1 << t[1][1] - t[2][1], t[2][1] - t[0][1], t[0][1] - t[1][1],
//		t[2][0] - t[1][0], t[0][0] - t[2][0], t[1][0] - t[0][0];
//	B1 = B1 / (2 * t.getArea());
//	Eigen::Matrix<Real,3,6> B2;
//	B2 << 4*L[0]-1, 0       , 0       , 0        , 4*L[2], 4*L[1],
//		  0       , 4*L[1]-1, 0       , 4*L[2]   , 0     , 4*L[0],
//		  0       , 0       , 4*L[2]-1, 4*L[1]   , 4*L[0], 0     ;
//	return(B1*B2*coefficients);
//}

//template Real evaluate_point<1>(const Triangle<3>& t, const Point& point, const Eigen::Matrix<Real,3,1>& coefficients);
//template Real evaluate_point<2>(const Triangle<6>& t, const Point& point, const Eigen::Matrix<Real,6,1>& coefficients);
//
//template Eigen::Matrix<Real,2,1> evaluate_der_point<1>(const Triangle<3>& t, const Point& point, const Eigen::Matrix<Real,3,1>& coefficients);
//template Eigen::Matrix<Real,2,1> evaluate_der_point<2>(const Triangle<6>& t, const Point& point, const Eigen::Matrix<Real,6,1>& coefficients);

