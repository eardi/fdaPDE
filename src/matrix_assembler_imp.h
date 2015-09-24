#ifndef MATRIX_ASSEMBLER_IMP_H_
#define MATRIX_ASSEMBLER_IMP_H_


//template<UInt ORDER, typename Integrator, typename P, typename A>
//void Assembler::operKernel(EOExpr<P, A> oper,const MeshHandler<ORDER>& mesh,
//	                     FiniteElement<Integrator, ORDER>& fe)
//	                     {
//  	UInt ntria=mesh.num_triangles(),
//		 nelem=mesh.num_nodes();
//
//  	std::vector<coeff> coefflist(ntria*(Integrator::NNODES)*(Integrator::NNODES)); // total number of summing terms
//  	//oper_mat_.resize(nelem,nelem);
//  	SPoper_mat_.resize(nelem,nelem);
//
//  	for(auto i=0; i<mesh.num_triangles(); i++){
//
//	    fe.updateElement(mesh.getTriangle(i));
//
//		// Vector of vertices indices
//		Eigen::Matrix<UInt,ORDER*3,1> identifiers;
//
//		//create a vector of the total number of contributes
//		for( auto q=0; q<ORDER*3; q++)
//		identifiers(q)=mesh.getTriangle(i)[q].id();
//
//		//localM=localMassMatrix(currentelem);
//		for(int i = 0; i < 3*ORDER; i++)
//		{
//			for(int j = 0; j < 3*ORDER; j++)
//			{
//				Real s=0;
//
//				for(int l = 0;l < Integrator::NNODES; l++)
//				{
//					s += oper(fe,i,j,l) * fe.getDet() * fe.getAreaReference()* Integrator::WEIGHTS[l];//(*)
//				}
//			  coefflist.emplace_back(coeff(identifiers(i),identifiers(j),s));
//			  //oper_mat_(identifiers(j),identifiers(k))+=oper(j,k);
//			}
//		}
//
//		}
//
//	SPoper_mat_.setFromTriplets(coefflist.begin(),coefflist.end());
//	//cout<<"done!"<<endl;;
//}

template<UInt ORDER, typename Integrator, typename A>
void Assembler::operKernel(EOExpr<A> oper,const MeshHandler<ORDER>& mesh,
	                     FiniteElement<Integrator, ORDER>& fe, SpMat& OpMat)
{
	std::vector<coeff> triplets;


  	for(auto t=0; t<mesh.num_triangles(); t++)
  	{
		fe.updateElement(mesh.getTriangle(t));

		// Vector of vertices indices (link local to global indexing system)
		std::vector<UInt> identifiers;
		identifiers.resize(ORDER*3);
		for( auto q=0; q<ORDER*3; q++)
			identifiers[q]=mesh.getTriangle(t)[q].id();


		//localM=localMassMatrix(currentelem);
		for(int i = 0; i < 3*ORDER; i++)
		{
			for(int j = 0; j < 3*ORDER; j++)
			{
				Real s=0;

				for(int l = 0;l < Integrator::NNODES; l++)
				{
					s += oper(fe,i,j,l) * fe.getDet() * fe.getAreaReference()* Integrator::WEIGHTS[l];//(*)
					//std::cout<<"("<<i<<","<<j<<","<<l<<"): "<<oper(fe,i,j,l)<< " " <<fe.getDet() << " " << fe.getAreaReference()<< " " << Integrator::WEIGHTS[l]<<"\n";
				}
			  triplets.push_back(coeff(identifiers[i],identifiers[j],s));
			}
		}

	}

  	UInt nnodes = mesh.num_nodes();
  	OpMat.resize(nnodes, nnodes);
	OpMat.setFromTriplets(triplets.begin(),triplets.end());
	//cout<<"done!"<<endl;;
}

template<UInt ORDER, typename Integrator>
void Assembler::forcingTerm(const MeshHandler<ORDER>& mesh,
	                     FiniteElement<Integrator, ORDER>& fe, const ForcingTerm& u, VectorXr& forcingTerm)
{

	forcingTerm = VectorXr::Zero(mesh.num_nodes());

  	for(auto t=0; t<mesh.num_triangles(); t++)
  	{
		fe.updateElement(mesh.getTriangle(t));

		// Vector of vertices indices (link local to global indexing system)
		std::vector<UInt> identifiers;
				identifiers.resize(ORDER*3);

		for( auto q=0; q<ORDER*3; q++)
			identifiers[q]=mesh.getTriangle(t)[q].id();


		//localM=localMassMatrix(currentelem);
		for(int i = 0; i < 3*ORDER; i++)
		{
			Real s=0;

			for(int iq = 0;iq < Integrator::NNODES; iq++)
			{
				UInt globalIndex = fe.getGlobalIndex(iq);
				s +=  fe.phiMaster(i,iq)* u(globalIndex) * fe.getDet() * fe.getAreaReference()* Integrator::WEIGHTS[iq];//(*)
			}
			forcingTerm[identifiers[i]] += s;
		}

	}
	//cout<<"done!"<<endl;;
}
    
    
#endif
