
#include <iostream>
#include "../../src/fdaPDE.h"
//#include "IO_handler.hpp"
#include "../../src/mesh_objects.h"
#include "../../src/mesh.h"
#include "../../src/finite_element.h"
#include "../../src/integration.h"
#include "../../src/matrix_assembler.h"
#include "../../src/param_functors.h"
#include "../../src/solver.h"

#include "../../src/regressionData.h"
#include "../../src/mixedFERegression.h"

//#include <iomanip>

//#include "regression_PDE.hpp"
#include "../../src/evaluator.h"
//
int main()
{
	//Simple Mesh
	//Real points[] = {0, 0, 0.5, 1, 1, 0, 1, 0.5, 1, 0};
	Real points[] = {0, 0, 0.5, 1, 1, 0, 0.25, 0.25, 0.75, 1, 0.75, 0.5, 0, 1, 0.5, 1, 0, 0.5,0.25, 0.75, 0.25, 0.5, 0.75, 0};
	UInt edges[] = {1,0,2,2,4,3,0,0,2,1,4,3,2,4};
	UInt triangles[] = {1,2,4,0,4,2,2,3,0,6,9,6,7,10,11,5,8,8};
	UInt neighbors[] = {2,-1,0,-1,-1,-1,-1,2,1};
	UInt num_nodes = 12;
	UInt num_edges = 7;
	UInt num_triangles = 3;

	MeshHandler<2> mesh(points, edges, triangles, neighbors, num_nodes, num_edges, num_triangles);


	//Simple regression problem
	std::vector<Point> locations;

	VectorXr observations(5);
	observations << 1,2,1,2,1;

	//MatrixXr covariates;
	MatrixXr covariates(5,1);
	covariates << 1, 2, 3, 4, 5;

	UInt order = 1;
	std::vector<Real> lambda(5);
	lambda[0] = 1; lambda[1] = 2; lambda[2] = 3; lambda[3] = 4; lambda[4] = 5;

	Eigen::Matrix<Real,2,2> K0;
		K0 << 1,0,0,1;
	std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > >
		K (3*IntegratorTriangleP4::NNODES, K0);

	Eigen::Matrix<Real,2,1> beta0;
		beta0 << 0,0;
	std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > >
		beta (3*IntegratorTriangleP4::NNODES, beta0);

	std::vector<Real> c (3*IntegratorTriangleP4::NNODES,0);
	std::vector<Real> u (3*IntegratorTriangleP4::NNODES,0);

	std::vector<UInt> dirichlet_indices;
	std::vector<Real> dirichlet_values;

	RegressionDataEllipticSpaceVarying regressionData(locations, observations, order, lambda, K, beta, c, u, covariates, dirichlet_indices, dirichlet_values, false);

	MixedFERegression<RegressionDataEllipticSpaceVarying, IntegratorTriangleP4, 2> regression(mesh,regressionData);
	regression.smoothEllipticPDESpaceVarying();

	std::cout<<"Solution"<<std::endl;
	for(auto i = 0; i< regression.getSolution().size();++i)
	{
		std::cout<<regression.getSolution()[i];
	}

	Real result[2];
	bool isinside[2];
	Real X[1]; X[0]=-10;
	Real Y[1]; Y[0]=10;
	UInt n_X = 1;
	Evaluator<2> evaluator(mesh);
	const Real* coeff = regression.getSolution()[0].data();

	evaluator.eval(X, Y, n_X, coeff, order, false, result, isinside);
	std::cout<<"\n Evaluation "<<result[0];
	std::cout<<"\n Is Inside "<<isinside[0];

	return 0;
}


//int main()
//{
//	//Simple Mesh
//	//Real points[] = {0, 0, 0.5, 1, 1, 0, 1, 0.5, 1, 0};
//	Real points[] = {0, 0, 0.5, 1, 1, 0, 0.25, 0.25, 0.75, 1, 0.75, 0.5, 0, 1, 0.5, 1, 0, 0.5,0.25, 0.75, 0.25, 0.5, 0.75, 0};
//	UInt edges[] = {1,0,2,2,4,3,0,0,2,1,4,3,2,4};
//	UInt triangles[] = {1,2,4,0,4,2,2,3,0,6,9,6,7,10,11,5,8,8};
//	UInt neighbors[] = {2,-1,0,-1,-1,-1,-1,2,1};
//	UInt num_nodes = 12;
//	UInt num_edges = 7;
//	UInt num_triangles = 3;
//
//	MeshHandler<2> mesh(points, edges, triangles, neighbors, num_nodes, num_edges, num_triangles);
//
//	//Simple regression problem
//	std::vector<Point> locations(5);
//	locations[0]= Point(0,0); locations[1] = Point(0,1); locations[2] = Point(0.5,0.5); locations[3] = Point(1,1); locations[4] = Point(1,0);
//
//	VectorXr observations(5);
//	observations << 1,2,1,2,1;
//
//	//MatrixXr covariates;
//	MatrixXr covariates(5,1);
//	covariates << 1, 2, 3, 4, 5;
//
//		//std::cout<<covariates<<std::endl;
//
//	UInt order = 1;
//	std::vector<Real> lambda(5);
//	lambda[0] = 1; lambda[1] = 2; lambda[2] = 3; lambda[3] = 4; lambda[4] = 5;
//	Eigen::Matrix<Real,2,2> K;
//		K << 1,0,0,1;
//		std::cout<< K ;
//	Eigen::Matrix<Real,2,1> beta;
//		beta << 0,0;
//	Real c = 0;
//	std::vector<UInt> dirichlet_indices;
//	std::vector<Real> dirichlet_values;
//
//	RegressionDataElliptic regressionData(locations, observations, order, lambda, K, beta, c, covariates, dirichlet_indices, dirichlet_values, true);
//
//	MixedFERegression<RegressionDataElliptic, IntegratorTriangleP4, 2> regression(mesh,regressionData);
//	regression.smoothEllipticPDE();
//
//	std::cout<<"Solution"<<std::endl;
//	for(auto i = 0; i< regression.getSolution().size();++i)
//	{
//		std::cout<<regression.getSolution()[i];
//		std::cout<<"DOF "<<regression.getDOF()[i];
//	}
//	return 0;
//}


//int main()
//{
//	//Simple Mesh
//	//Real points[] = {0, 0, 0.5, 1, 1, 0, 1, 0.5, 1, 0};
//	Real points[] = {0, 0, 0.5, 1, 1, 0, 0.25, 0.25, 0.75, 1, 0.75, 0.5, 0, 1, 0.5, 1, 0, 0.5,0.25, 0.75, 0.25, 0.5, 0.75, 0};
//	UInt edges[] = {1,0,2,2,4,3,0,0,2,1,4,3,2,4};
//	UInt triangles[] = {1,2,4,0,4,2,2,3,0,6,9,6,7,10,11,5,8,8};
//	UInt neighbors[] = {2,-1,0,-1,-1,-1,-1,2,1};
//	UInt num_nodes = 12;
//	//UInt num_nodes = 5
//	UInt num_edges = 7;
//	UInt num_triangles = 3;
//
//	MeshHandler<2> mesh(points, edges, triangles, neighbors, num_nodes, num_edges, num_triangles);
//
//	//Simple regression problem
//	std::vector<Point> locations;
//	//locations[0]= Point(0,0); locations[1] = Point(0,1); locations[2] = Point(0.5,0.5); locations[3] = Point(1,1); locations[4] = Point(1,0);
//
//	VectorXr observations(5);
//	observations << 1,2,1,2,1;
//
//	//MatrixXr covariates;
//	MatrixXr covariates(5,1);
//	covariates << 1, 2, 3, 4, 5;
//
//		//std::cout<<covariates<<std::endl;
//
//	UInt order = 2;
//	std::vector<Real> lambda(3);
//	lambda[0] = 1; lambda[1] = 2; lambda[2] = 3;
//	std::vector<UInt> dirichlet_indices;
//	std::vector<Real> dirichlet_values;
//
//	RegressionData regressionData(locations, observations, order, lambda, covariates, dirichlet_indices, dirichlet_values);
//
//
//	MixedFERegression<RegressionData, IntegratorTriangleP4, 2> regression(mesh,regressionData);
//	regression.smoothLaplace();
//
//	std::cout<<"Solution"<<std::endl;
//	for(auto i = 0; i< regression.getSolution().size();++i)
//	{
//		std::cout<<regression.getSolution()[i];
//	}
//
//	return 0;
//}


//int main()
//{
//	VectorXr observations(5);
//	observations << 1,2,1,2,1;
//	MatrixXr design_matrix(5,1);
//	design_matrix << 1, 2, 3, 4, 5;
//	//std::cout<<design_matrix<<std::endl;
//	UInt order = 1;
//	std::vector<UInt> dirichlet_indexes;
//	std::vector<Real> dirichlet_values;
//
//	//Real points[] = {0, 0, 0.5, 1, 1, 0, 1, 0.5, 1, 0};
//	Real points[] = {0, 0, 0.5, 1, 1, 0, 0.25, 0.25, 0.75, 1, 0.75, 0.5, 0, 1, 0.5, 1, 0, 0.5,0.25, 0.75, 0.25, 0.5, 0.75, 0};
//	UInt edges[] = {1,0,2,2,4,3,0,0,2,1,4,3,2,4};
//	UInt triangles[] = {1,2,4,0,4,2,2,3,0,6,9,6,7,10,11,5,8,8};
//	UInt neighbors[] = {2,-1,0,-1,-1,-1,-1,2,1};
//	UInt num_nodes = 12;
//	UInt num_edges = 7;
//	UInt num_triangles = 3;
//
//	MeshHandler<2> mesh(points, edges, triangles, neighbors, num_nodes, num_edges, num_triangles);
//
//	//Simple regression problem
//	std::vector<Point> locations(5);
//	locations[0]= Point(0,0); locations[1] = Point(0,1); locations[2] = Point(0.5,0.5); locations[3] = Point(1,1); locations[4] = Point(1,0);
//
//	//MatrixXr covariates;
//	MatrixXr covariates(5,1);
//	covariates << 1, 2, 3, 4, 5;
//
//		//std::cout<<covariates<<std::endl;
//
//	std::vector<Real> lambda(1);
//	lambda[0] = 1;
//	std::vector<UInt> dirichlet_indices;
//
//	Eigen::Matrix<Real,2,2> K;
//	K << 2,0,0,2;
//	Eigen::Matrix<Real,2,1> beta;
//	beta << 0,0;
//	Real mu = 1;
//	Eigen::Matrix<Real,2,1> d;
//	d << 1,1;
//	Real c = 0;
//
//	RegressionDataElliptic regressionData(locations, observations, order, lambda, K, beta, c, covariates, dirichlet_indices, dirichlet_values, true);
//
//
//	FiniteElement<IntegratorTriangleP4,2> fe;
//
//	Diffusivity diff;
//
//	typedef EOExpr<Mass> ETMass;
//	typedef EOExpr<Stiff> ETStiff;
//	typedef EOExpr<Grad> ETGrad;
//
//	Mass EMass;
//	Stiff EStiff;
//	Grad EGrad;
//
//	ETMass mass(EMass);
//	ETStiff stiff(EStiff);
//	ETGrad grad(EGrad);
//
//	SpMat AMat;
//
//	Assembler::operKernel(stiff[K],mesh,fe, AMat);
//
//	std::cout<<AMat;

//	Model<IntegratorTriangleP4, 2> model(mesh,fe);
//
//	VectorXr forcing_term(mesh.num_nodes());
//	VectorXr rightside(12);
//	rightside << 1,2,1,2,1,2,1,2,1,2,1,2;
//	std::vector<coeff> tripletsData(12);
//	for (int i =0; i<12;++i)
//		tripletsData[i] = coeff(i,i,1);
//
////(-iohandler.getLambda())
//	model.build(tripletsData,  (-iohandler.getLambda())*stiff,  mu*mass, rightside, forcing_term);
//	model.solve<SpLU>();
//
//	std::cout<<model.getSolution();
//
//	return 0;
//}


//int main()
//{
//	VectorXr observations(5);
//	observations << 1,2,1,2,1;
//	MatrixXr design_matrix(5,1);
//	design_matrix << 1, 2, 3, 4, 5;
//	//std::cout<<design_matrix<<std::endl;
//	UInt order = 1;
//	Real lambda = 1;
//	std::vector<UInt> dirichlet_indexes;
//	std::vector<Real> dirichlet_values;
//
//	//Real points[] = {0, 0, 0.5, 1, 1, 0, 1, 0.5, 1, 0};
//	Real points[] = {0, 0, 0.5, 1, 1, 0, 0.25, 0.25, 0.75, 1, 0.75, 0.5, 0, 1, 0.5, 1, 0, 0.5,0.25, 0.75, 0.25, 0.5, 0.75, 0};
//	UInt edges[] = {1,0,2,2,4,3,0,0,2,1,4,3,2,4};
//	UInt triangles[] = {1,2,4,0,4,2,2,3,0,6,9,6,7,10,11,5,8,8};
//	UInt neighbors[] = {2,-1,0,-1,-1,-1,-1,2,1};
//	UInt num_nodes = 12;
//	UInt num_edges = 7;
//	UInt num_triangles = 3;
//
//	MeshHandler<2> mesh(points, edges, triangles, neighbors, num_nodes, num_edges, num_triangles);
//	//IOHandler iohandler(observations, design_matrix, order, lambda, dirichlet_indexes, dirichlet_values);
//
//	Triangle<3> t3 = Triangle<3>(Id(1), std::vector<Point> { {Point(0,0,0,1), Point(1,0,0,0),Point(2,0,0.5,0.5)} });
//	Triangle<6> t6 = Triangle<6>(Id(2), std::vector<Point> { {Point(0,0,0,1), Point(1,0,0,0),Point(2,0,0.5,0.5), Point(3,0,0.25,0.25), Point(4,0,0.25,0.75),Point(5,0,0,0.5)} });
//
//
//
//	t3.print(std::cout);
//	std::cout<<t3.getArea();
//	t6.print(std::cout);
//
//
//	Eigen::Matrix<Real,3,1> coefficientsa0;
//	Eigen::Matrix<Real,3,1> coefficientsa1;
//	Eigen::Matrix<Real,3,1> coefficientsa2;
//
//	Eigen::Matrix<Real,6,1> coefficientsb0;
//	Eigen::Matrix<Real,6,1> coefficientsb1;
//	Eigen::Matrix<Real,6,1> coefficientsb2;
//	Eigen::Matrix<Real,6,1> coefficientsb3;
//	Eigen::Matrix<Real,6,1> coefficientsb4;
//	Eigen::Matrix<Real,6,1> coefficientsb5;
//
//	coefficientsa0<< 1, 0, 0;
//	coefficientsa1<< 0, 1, 0;
//	coefficientsa2<< 0, 0, 1;
//
//	coefficientsb0<< 1, 0, 0, 0, 0, 0;
//	coefficientsb1<< 0, 1, 0, 0, 0, 0;
//	coefficientsb2<< 0, 0, 1, 0, 0, 0;
//	coefficientsb3<< 0, 0, 0, 1, 0, 0;
//	coefficientsb4<< 0, 0, 0, 0, 1, 0;
//	coefficientsb5<< 0, 0, 0, 0, 0, 1;
//
//	Real res = evaluate_point<1>(t3,Point(0.25,0.25),coefficientsa0);
//
//	std::cout<<"Eval 1: "<<evaluate_point<1>(t3,Point(0.25,0.25),coefficientsa0)<<std::endl;
//
//	std::cout<<"Eval 1 der: "<<evaluate_der_point<1>(t3,Point(0.25,0.25),coefficientsa0)<<std::endl;
//	std::cout<<"Eval 2 der: "<<evaluate_der_point<1>(t3,Point(0.25,0.25),coefficientsa1)<<std::endl;
//	std::cout<<"Eval 3 der: "<<evaluate_der_point<1>(t3,Point(0.25,0.25),coefficientsa2)<<std::endl;
//
//	std::cout<<"Eval 1: "<<evaluate_point<2>(t6,Point(0.25,0.25),coefficientsb0)<<std::endl;
//
//	std::cout<<"Eval 1 der: "<<evaluate_der_point<2>(t6,Point(0.25,0.25),coefficientsb0)<<std::endl;
//	std::cout<<"Eval 2 der: "<<evaluate_der_point<2>(t6,Point(0.25,0.25),coefficientsb1)<<std::endl;
//	std::cout<<"Eval 3 der: "<<evaluate_der_point<2>(t6,Point(0.25,0.25),coefficientsb2)<<std::endl;
//	std::cout<<"Eval 4 der: "<<evaluate_der_point<2>(t6,Point(0.25,0.25),coefficientsb3)<<std::endl;
//	std::cout<<"Eval 5 der: "<<evaluate_der_point<2>(t6,Point(0.25,0.25),coefficientsb4)<<std::endl;
//	std::cout<<"Eval 6 der: "<<evaluate_der_point<2>(t6,Point(0.25,0.25),coefficientsb5)<<std::endl;
//
//
//	FiniteElement<IntegratorTriangleP4,2> fe;
//	fe.updateElement(t6);
//
//	Point p = fe.coorQuadPt(0);
//	std::cout<<"PROVA TRASF: ("<<p[0]<<","<<p[1]<<")"<<std::endl;
//
//	std::cout<<"J^-1 ref "<<std::endl<<t3.getM_invJ()<<std::endl<<"Cross_Prod"<<std::endl;
//	for (int i =0; i<6; i++)
//	{
//		for(int j=0; j<6; j++)
//		{
//			Real s=0;
//			for (int iq=0;iq<IntegratorTriangleP4::NNODES; iq++)
//				s+=fe.phiMaster(i,iq)*fe.phiMaster(j,iq)*IntegratorTriangleP4::WEIGHTS[iq]*(1./2)*fe.getDet();
//		std::cout<<s<<"\t";
//		}
//		std::cout<<std::endl;
//	}
//
//	fe.updateElement(t6);
//	for (int i =0; i<6; i++)
//		{
//			for(int j=0; j<6; j++)
//			{
//				Real s1=0,s2=0;
//				for (int iq=0;iq<3; iq++)
//				{
//					s1 += fe.invTrJPhiDerMaster(i,0,iq)*fe.invTrJPhiDerMaster(j,0,iq)+
//							fe.invTrJPhiDerMaster(i,1,iq)*fe.invTrJPhiDerMaster(j,1,iq);
//				}
//			std::cout<<"("<<s1<<","<<s2<<")"<<"\t";
//			}
//		std::cout<<std::endl;
//	}
//
//	//Diffusivity diff;
//
//	typedef EOExpr<Mass> ETMass;
//	typedef EOExpr<Stiff> ETStiff;
//	typedef EOExpr<Grad> ETGrad;
//
//	Mass EMass;
//	Stiff EStiff;
//	Grad EGrad;
//
//	ETMass mass(EMass);
//	ETStiff stiff(EStiff);
//	ETGrad grad(EGrad);
//
//	Assembler assembleMass;
//
//	Eigen::Matrix<Real,2,2> K;
//	K << 2,0,0,2;
//	Real mu = 1;
//	Eigen::Matrix<Real,2,1> d;
//	d << 1,1;
//
//	SpMat Mass;
//	assembleMass.operKernel(mu*mass + stiff[K] + dot(d,grad),mesh,fe, Mass);
//
//	std::cout<<"MASS"<<std::endl;
//	std::cout<< Mass;
//
//	//RegressionPDE<ExactFirstOrder,1> regression(mesh,iohandler);
//	//const VectorXr& result_eigen = regression.smoothBase();
//
//	//std::cout<<result_eigen;
//
//	return 0;
//}
