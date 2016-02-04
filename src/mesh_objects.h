#ifndef __MESH_OBJECTS_HPP__
#define __MESH_OBJECTS_HPP__


#include "fdaPDE.h"

//Accord the NotValid meaning value
//const UInt NVAL=std::numeric_limits<UInt>::max();

typedef UInt Id;
typedef UInt BcId;

//!  This class gives some common methods to all mesh objects.
class Identifier{
public:
	
	//! An static const Unisgned Integer.
    /*! Needed to identify the Not Valid Id. */
	static const UInt NVAL;
	//Identifier():id_(NVAL),bcId_(NVAL){}
	Identifier(UInt id):id_(id),bcId_(NVAL){}
	Identifier(UInt id, UInt bcId):id_(id),bcId_(bcId){}
	
	bool unassignedId()const {return id_==NVAL;}
	bool unassignedBc()const {return bcId_==NVAL;}
	
	Id id() const {return id_;}
	BcId bcId() const {return bcId_;}
	Id getId() const {return id_;}
	
	
	protected:
	Id id_;
	BcId bcId_;
};


//!  This class implements a 2D point
class Point: public Identifier{
public:
	
	static const UInt ndim = 2;

	Point(): Identifier(NVAL, NVAL){coord_.resize(2);};
	Point(Real x, Real y):Identifier(NVAL, NVAL)
		{coord_.resize(2);coord_[0]=x; coord_[1]=y;}
	Point(Id id, BcId bcId, Real x, Real y):Identifier(id, bcId)
		{coord_.resize(2);coord_[0]=x; coord_[1]=y;}
	void print(std::ostream & out) const;
	Real operator[](UInt i) const {return coord_[i];}
private:
	std::vector<Real> coord_;
	//std::array<Real, 2> coord_;
};


//!  This class implements an Edge, as an objects composed by two 2D points.
class Edge: public Identifier{
  public:
    static const UInt NNODES=2;
    static const UInt numSides=1;
    static const UInt myDim=1;
    
    Edge():Identifier(NVAL, NVAL){points_.resize(2);};
    Edge(Id id, BcId bcId, const Point& start,const Point& end):Identifier(id, bcId)
    {points_.resize(2); points_[0] = start; points_[1] = end;}
    
    void print(std::ostream & out) const;
    Point getFirst() const {return points_[0];}
    Point getEnd() const {return points_[1];}
    
    
    Point operator[](UInt i) const {return points_[i];}
   
 private:
	// I don't store directly a eigen matrix because of the limitations
	// of the current problems of alignement (see eigen 3.0 documentation)
	// It is not very efficient and should be changed asap
    //std::array<Point, NNODES> points_;
    std::vector<Point> points_;
    //std::array<std::reference_wrapper<Point>, NNODES> M_points;
  };
  
//!  This class implements a Triangle as an objects composed by three or six nodes.
/*!
 *  The first three nodes represent the vertices, the others the internal nodes,
 *  following this enumeration: !IMPORTANT! different from Sangalli code!
 * 
 * 				3
 * 				*
 * 			 /	   \
 * 		  5 *		* 4 
 * 		  /			  \
 * 		 *______*______*		
 * 		1		6		2
*/

template <UInt NNODES>
class Triangle : public Identifier {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    static const UInt numVertices=3;
    static const UInt numSides=3;
	static const UInt myDim=2;
    
    //! This constructor creates an "empty" Triangle, with an Id Not Valid
	Triangle():Identifier(NVAL){points_.resize(NNODES);}
	
	//! This constructor creates a Triangle, given its Id and an std array with the three object Point the will define the Triangle
    Triangle(Id id, const std::vector<Point> points) : Identifier(id)
	{ this->init(points); }
	
	//! Overloading of the operator [],  taking the Node number and returning a node as Point object.
    /*!
     * For node numbering convention see: 
      \param i an integer argument.
      \return the Point object
    */
	Point operator[](UInt i) const {return points_[i];}
	
	//! A member that computes the barycentric coordinates.
    /*!
      \param point a Point object
      \return The three baricentric coordinates of the point
    */

	Real getDetJ() const {return detJ_;}
	Eigen::Matrix<Real,2,2> getM_J() const {return M_J_;}
	Eigen::Matrix<Real,2,2> getM_invJ() const {return M_invJ_;}
	Eigen::Matrix<Real,2,2> getMetric() const {return metric_;}
	//! A member returning the area of the finite element
	    /*!
	      \return a Real value representing the area of the triangle from which we updated the element
	      \sa  updateElement(Triangle<Integrator::NNODES> t)
	    */
	Real getArea() const {return (0.5 * detJ_);}

	Eigen::Matrix<Real,3,1> getBaryCoordinates(const Point& point) const;
	
	//! A member that tests if a Point is located inside a Triangle.
    /*!
      \param point a Point object.
      \return True if the point is inside the triangle
    */
	bool isPointInside(const Point& point) const;
	
	//! A memeber that verifies which edge separates the Triangle from a Point.
    /*!
      \param point a Point object.
      \return The number of the Edge that separates the point
      from the triangle and -1 if the point is inside the triangle.
    */
	int getPointDirection(const Point& point) const;
	
	//! A member that prints the main properties of the triangle
    /*!
      \param out a std::outstream.
    */
	void print(std::ostream & out) const;

private:
	//std::array<Point, NNODES> points_;
	std::vector<Point> points_;
	Eigen::Matrix<Real,2,2> M_J_;
	Eigen::Matrix<Real,2,2> M_invJ_;
	Eigen::Matrix<Real,2,2> metric_;
	Real detJ_;
	void init(const std::vector<Point> &points);
};

template <UInt NNODES>
const int Triangle<NNODES>::myDim;

//! A function for the evaluation of point value in a triangle.
/*!
  \param t a Triangle object
  \param point a point object
  \param coefficients a Eigen vector specifing the coefficients of the Lagrangian
		 base (1st or 2nd order) defined on the Triangle.
  \return The point evaluation of the function defined by the coefficients on
  the triangle
    */

//template <UInt ORDER>
//Eigen::Matrix<Real,2,1> evaluate_der_point(const Triangle<3*ORDER>& t, const Point& point, const Eigen::Matrix<Real,3*ORDER,1>& coefficients)
//{
//	//std::cerr<< "TRYING TO EVALUATE ORDER NOT IMPLEMENTED" << std::endl;
//	Eigen::Matrix<Real,2,1> null;
//	return(null);
//}

//template <UInt ORDER>
//Real evaluate_point(const Triangle<3*ORDER>& t, const Point& point, const Eigen::Matrix<Real,3*ORDER,1>& coefficients);

//template <>
//Real evaluate_point<1>(const Triangle<3>& t, const Point& point, const Eigen::Matrix<Real,3,1>& coefficients);
//
//template <>
//Real evaluate_point<2>(const Triangle<6>& t, const Point& point, const Eigen::Matrix<Real,6,1>& coefficients);

//template <UInt ORDER>
//Eigen::Matrix<Real,2,1> evaluate_der_point(const Triangle<3*ORDER>& t, const Point& point, const Eigen::Matrix<Real,3*ORDER,1>& coefficients);
//
//template <>
//Eigen::Matrix<Real,2,1> evaluate_der_point<1>(const Triangle<3>& t, const Point& point, const Eigen::Matrix<Real,3,1>& coefficients);

//template <>
//Eigen::Matrix<Real,2,1> evaluate_der_point<2>(const Triangle<6>& t, const Point& point, const Eigen::Matrix<Real,6,1>& coefficients);


template <UInt ORDER>
inline Real evaluate_point(const Triangle<3*ORDER>& t, const Point& point, const Eigen::Matrix<Real,3*ORDER,1>& coefficients)
{
	//std::cerr<< "TRYING TO EVALUATE ORDER NOT IMPLEMENTED" << std::endl;
	return 0;
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

template <UInt ORDER>
inline Eigen::Matrix<Real,2,1> evaluate_der_point(const Triangle<3*ORDER>& t, const Point& point, const Eigen::Matrix<Real,3*ORDER,1>& coefficients)
{
	//std::cerr<< "TRYING TO EVALUATE ORDER NOT IMPLEMENTED" << std::endl;
	Eigen::Matrix<Real,2,1> null;
	return(null);
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

#include "mesh_objects_imp.h"
#endif
