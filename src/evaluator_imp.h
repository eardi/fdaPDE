#ifndef __EVALUATOR_IMP_HPP__
#define __EVALUATOR_IMP_HPP__

template <UInt ORDER>
void Evaluator<ORDER>::eval(Real* X, Real *Y, UInt length, const Real *coef, UInt order, bool fast, Real* result, std::vector<bool>& isinside)
{
	
	
	Triangle<3*ORDER> current_triangle;
	std::vector<Triangle<3*ORDER> > starting_triangles;
	starting_triangles.resize(1);

	Point current_point;
	Eigen::Matrix<Real,3*ORDER,1> coefficients;
	
	starting_triangles[0] = mesh_.getTriangle(0);
	for (int i = 0; i<length; ++i)
	{
		current_point = Point(X[i],Y[i]);
		//current_triangle = mesh_.findLocationNaive(current_point);
		//std::cout<<"Looking For Position Walking... \n";
		current_triangle = mesh_.findLocationWalking(current_point, starting_triangles);
		//current_triangle.print(cout);
		//cout<<"triangle: "<< current_triangle.getId()<<endl;
		if(current_triangle.getId() == Identifier::NVAL && fast == false)
		{
			//To avoid problems with non convex mesh
			//std::cout<<"Position Not Found Walking... \n";
			current_triangle = mesh_.findLocationNaive(current_point);
		}
		if(current_triangle.getId() == Identifier::NVAL)         
		{
			//std::cout<<"Position Not Found Naively... \n";
			isinside[i]=false;
		}
		else 
		{
			isinside[i]=true;
			for (int i=0; i<(3*ORDER); ++i)		
			{
				coefficients[i] = coef[current_triangle[i].getId()];
			}
			result[i] = evaluate_point<ORDER>(current_triangle, current_point, coefficients);
			starting_triangles[0] = current_triangle;
		}
	}
}


#endif
