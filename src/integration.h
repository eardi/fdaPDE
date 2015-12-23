#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__

#include "fdaPDE.h"
#include "mesh_objects.h"



class IntegratorTriangleP2{
	public:
	static const UInt ORDER = 1;
	//Number of nodes
	static const UInt NNODES = 3;
	//Point locations
	static const std::vector<Point> NODES;
	static const std::vector<Real> WEIGHTS;
	};

class IntegratorTriangleP4{
	public:
	static const UInt ORDER = 2;
	//Number of nodes
	static const UInt NNODES = 6;
	//Point locations
	static const std::vector<Point> NODES;
	static const std::vector<Real> WEIGHTS;
	};

//#include "integration_imp.hpp"
#endif

