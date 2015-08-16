#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__

#include "FEMr.hpp"
#include "mesh_objects.hpp"



class IntegratorTriangleP2{
	public:
	static const UInt ORDER = 1;
	//Number of nodes
	static const UInt NNODES = 3;
	//Point locations
	static const std::vector<Point> NODES;
	static const std::vector<Real> WEIGHTS;
	};

const std::vector<Real> IntegratorTriangleP2::WEIGHTS = std::vector<Real>{ {1./3, 1./3, 1./3} };
const std::vector<Point> IntegratorTriangleP2::NODES = std::vector<Point> { {Point(1./6,1./6),Point(2./3,1./6),Point(1./6,2./3)} };

class IntegratorTriangleP4{
	public:
	static const UInt ORDER = 1;
	//Number of nodes
	static const UInt NNODES = 6;
	//Point locations
	static const std::vector<Point> NODES;
	static const std::vector<Real> WEIGHTS;
	};

const std::vector<Real> IntegratorTriangleP4::WEIGHTS = std::vector<Real>{ {
0.2233815897,
0.2233815897,
0.2233815897,
0.1099517437,
0.1099517437,
0.1099517437
} };

const std::vector<Point> IntegratorTriangleP4::NODES = std::vector<Point> { {
Point(0.4459484909,0.4459484909),
Point(0.4459484909,0.1081030182),
Point(0.1081030182,0.4459484909),
Point(0.0915762135, 0.0915762135),
Point(0.0915762135, 0.816847573),
Point(0.816847573, 0.0915762135)
} };

#endif
