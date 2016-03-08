#ifndef __EVALUATOR_IMP_HPP__
#define __EVALUATOR_IMP_HPP__

template <UInt ORDER>
void Evaluator<ORDER>::eval(Real* X, Real *Y, UInt length, const Real *coef, UInt order, bool redundancy, Real* result, std::vector<bool>& isinside)
{
	
	
	Triangle<3*ORDER> current_triangle;
	// std::vector<Triangle<3*ORDER> > starting_triangles; Problem with alignment not solved
	// by http://eigen.tuxfamily.org/dox-devel/group__TopicUnalignedArrayAssert.html
	// starting_triangles.resize(1);
	Triangle<3*ORDER> starting_triangle;


	Point current_point;
	Eigen::Matrix<Real,3*ORDER,1> coefficients;
	
	starting_triangle = mesh_.getTriangle(0);
	for (int i = 0; i<length; ++i)
	{
		current_point = Point(X[i],Y[i]);
		//current_triangle = mesh_.findLocationNaive(current_point);
		//std::cout<<"Finding point.. "<<i<<" from "<<current_triangle.getId()<<"\n";
		current_triangle = mesh_.findLocationWalking(current_point, starting_triangle);
		//current_triangle.print(cout);
		//std::cout<<"Walking...triangle: "<< current_triangle.getId()<<std::endl;
		if(current_triangle.getId() == Identifier::NVAL && redundancy == true)
		{
			//To avoid problems with non convex mesh
			//std::cout<<"Position Not Found Walking... \n";
			current_triangle = mesh_.findLocationNaive(current_point);
			//std::cout<<"Naively...triangle: "<< current_triangle.getId()<<std::endl;
		}
		if(current_triangle.getId() == Identifier::NVAL)         
		{
			//std::cout<<"Position Not Found Naively... \n";
			isinside[i]=false;
		}
		else 
		{
			isinside[i]=true;
			for (int j=0; j<(3*ORDER); ++j)
			{
				coefficients[j] = coef[current_triangle[j].getId()];
			}
			result[i] = evaluate_point<ORDER>(current_triangle, current_point, coefficients);
			starting_triangle = current_triangle;
		}
	}
}


#endif
