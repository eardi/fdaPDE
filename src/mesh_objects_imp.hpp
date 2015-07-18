#ifndef __MESH_OBJECTS_IMP_HPP__
#define __MESH_OBJECTS_IMP_HPP__

const UInt Identifier::NVAL=std::numeric_limits<UInt>::max();


template <UInt NNODES>
void Triangle<NNODES>::init(const std::vector<Point> &points)
{
	points_ = points;

	Triangle<NNODES> &t = *this;
	Point d1(t[1][0]-t[0][0], t[1][1]-t[0][1]);
	Point d2(t[2][0]-t[0][0], t[2][1]-t[0][1]);   //Point d2 = t[2] - t[0]; non funziona, reimplementare sottrazione

	M_J_(0,0) = d1[0];			// (x2-x1)
	M_J_(1,0) = d1[1];			// (y2-y1)
	M_J_(0,1) = d2[0];			// (x3-x1)
	M_J_(1,1) = d2[1];			// (y3-y1)

	detJ_ = M_J_(0,0) * M_J_(1,1) - M_J_(1,0) * M_J_(0,1);

	Real idet = 1. / detJ_;

	M_invJ_(0,0) =  idet * M_J_(1,1);	// (y3-y1)	-(x3-x1)
	M_invJ_(1,0) = -idet * M_J_(1,0);	// -(y2-y1) (x2-x1)
	M_invJ_(0,1) = -idet * M_J_(0,1);	//
	M_invJ_(1,1) =  idet * M_J_(0,0);	//	è la trasposta di quella della Sangalli (Ael)

	metric_ = M_invJ_ * M_invJ_.transpose();
}

template <UInt NNODES>
Eigen::Matrix<Real,3,1> Triangle<NNODES>::getBaryCoordinates(const Point& point) const{
	
	Triangle<NNODES> t=*this;
	Eigen::Matrix<Real,3,1> lambda;
	Eigen::Matrix<Real,4,1> bary_coef;
	//Real eps = 2.2204e-016,
	//	 tolerance = 10000 * eps;

	//cout<<"primovert  "<<t[0](0)<<endl;
	//cout<<"primovert  "<<t[0](1)<<endl;
		 
	bary_coef[0] = t[0][0]-t[2][0];  //x1-x3
	bary_coef[1] = t[1][0]-t[2][0];  //x2-x3
	bary_coef[2] = t[0][1]-t[2][1];  //y1-y3
	bary_coef[3] = t[1][1]-t[2][1];  //y2-y3
	//cout<<baryccoef<<endl;
	//cout<<baryccoef<<endl;
	
	Real detT = bary_coef[0]*bary_coef[3]-bary_coef[1]*bary_coef[2];
	bary_coef = bary_coef / detT;
	//cout<<"detT  "<<detT<<endl;
	
	//Compute barycentric coordinates for the point
	Real x_diff_third = point[0] - t[2][0],
		 y_diff_third = point[1] - t[2][1];
		 
	lambda[0] =  (bary_coef[3] * x_diff_third - bary_coef[1] * y_diff_third),
	lambda[1] = (-bary_coef[2] * x_diff_third + bary_coef[0] * y_diff_third),
	lambda[2] = 1 - lambda[0] - lambda[1];

	return lambda;

}  


template <UInt NNODES>
bool Triangle<NNODES>::isPointInside(const Point& point) const
{
	Real eps = 2.2204e-016,
		 tolerance = 10000 * eps;

	Eigen::Matrix<Real,3,1> lambda = getBaryCoordinates(point);
	
	return ((-tolerance <= lambda[0] && lambda[0] <= 1+tolerance) && 
			(-tolerance <= lambda[1] && lambda[1] <= 1+tolerance) &&
			(-tolerance <= lambda[2] && lambda[2] <= 1+tolerance) );
	
}  


// TO BE FIXED: if one dir -1, try with others
template <UInt NNODES>
int Triangle<NNODES>::getPointDirection(const Point& point) const
{
	Real eps = 2.2204e-016,
		 tolerance = 10000 * eps;

	Eigen::Matrix<Real,3,1> lambda = getBaryCoordinates(point);
	
	//Find the minimum coordinate (if negative stronger straight to the point searched)
	int min_index;
	if(lambda[0] < lambda[1] && lambda[0] < lambda[2]) 		min_index = 0;
	else if(lambda[1] < lambda[0] && lambda[1] < lambda[2]) min_index = 1;
	else 													min_index = 2;
	
	if(lambda[min_index] < -tolerance) 	return min_index;
	else 							   	return -1;
}

inline
void Point::print(std::ostream & out) const
{
	out<<"Point -"<< id_ <<"- "<<"("<<coord_[0]<<","<<coord_[1]<<")"<<std::endl<<"------"<<std::endl;	
}

inline
void Edge::print(std::ostream & out) const
{
	out<<"Edge -"<< id_ <<"- "<<"("<<points_[0].getId()<<","<<points_[1].getId()<<")"<<std::endl;	
}

template <UInt NNODES>
void Triangle<NNODES>::print(std::ostream & out) const
{
	out<<"Triangle -"<< id_ <<"- ";
	for (UInt i=0; i<NNODES; ++i)
		out<<points_[i].getId()<<"  ";
	out<<std::endl;	
}
 
template <UInt ORDER>
Real evaluate_point(const Triangle<3*ORDER>& t, const Point& point, const Eigen::Matrix<Real,3*ORDER,1>& coefficients)
{
	//std::cerr<< "TRYING TO EVALUATE ORDER NOT IMPLEMENTED" << std::endl;
	return 0;
}

template <UInt ORDER>
inline Eigen::Matrix<Real,2,1> evaluate_der_point(const Triangle<3*ORDER>& t, const Point& point, const Eigen::Matrix<Real,3*ORDER,1>& coefficients)
{
	//std::cerr<< "TRYING TO EVALUATE ORDER NOT IMPLEMENTED" << std::endl;
	Eigen::Matrix<Real,2,1> null;
	return(null);
}


template <>
inline Real evaluate_point<1>(const Triangle<3>& t, const Point& point, const Eigen::Matrix<Real,3,1>& coefficients)
{
	Eigen::Matrix<Real,3,1> bary_coeff = t.getBaryCoordinates(point);
	//std::cout<< "B-coord: "<<bary_coeff<<std::endl;
	
	return(coefficients.dot(bary_coeff));	
}

template <>
inline Real evaluate_point<2>(const Triangle<6>& t, const Point& point, const Eigen::Matrix<Real,6,1>& coefficients)
{
	Eigen::Matrix<Real,3,1> bary_coeff = t.getBaryCoordinates(point);
	return( coefficients[0]*(2*bary_coeff[0]*bary_coeff[0]- bary_coeff[0]) +
            coefficients[1]*(2*bary_coeff[1]*bary_coeff[1] - bary_coeff[1]) +
            coefficients[2]*(2*bary_coeff[2]*bary_coeff[2] - bary_coeff[2]) +
            coefficients[3]*(4*bary_coeff[1]* bary_coeff[2])    +
            coefficients[4]*(4*bary_coeff[2]* bary_coeff[0])    +
            coefficients[5]*(4*bary_coeff[0]* bary_coeff[1]) );
                           
}

template <>
inline Eigen::Matrix<Real,2,1> evaluate_der_point<1>(const Triangle<3>& t, const Point& point, const Eigen::Matrix<Real,3,1>& coefficients)
{
	Eigen::Matrix<Real,2,3> B1;
	B1 << t[1][1] - t[2][1], t[2][1] - t[0][1], t[0][1] - t[1][1],
		t[2][0] - t[1][0], t[0][0] - t[2][0], t[1][0] - t[0][0];
	B1 = B1 / (2 * t.getArea());

	return(B1*coefficients);

}

template <>
inline Eigen::Matrix<Real,2,1> evaluate_der_point<2>(const Triangle<6>& t, const Point& point, const Eigen::Matrix<Real,6,1>& coefficients)
{
	Eigen::Matrix<Real,3,1> L = t.getBaryCoordinates(point);
	Eigen::Matrix<Real,2,3> B1;
	B1 << t[1][1] - t[2][1], t[2][1] - t[0][1], t[0][1] - t[1][1],
		t[2][0] - t[1][0], t[0][0] - t[2][0], t[1][0] - t[0][0];
	B1 = B1 / (2 * t.getArea());
	Eigen::Matrix<Real,3,6> B2;
	B2 << 4*L[0]-1, 0       , 0       , 0        , 4*L[2], 4*L[1],
		  0       , 4*L[1]-1, 0       , 4*L[2]   , 0     , 4*L[0],
		  0       , 0       , 4*L[2]-1, 4*L[1]   , 4*L[0], 0     ;
	return(B1*B2*coefficients);
}

#endif
