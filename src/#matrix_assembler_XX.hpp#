#ifndef MATRIX_ASSEMBLER_HPP_
#define MATRIX_ASSEMBLER_HPP_


#include "mesh_objects.hpp"
#include "finite_element.hpp"
#include "param_functors.hpp"
#include "RPDE.hpp"
//! A Stiff class: a class for the stiffness operator.

template <class Integrator, UInt ORDER>
class StiffXX{
  private:

  public:
	inline Real operator() (FiniteElement<Integrator, ORDER>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
	{
	   	Real s = 0;
	   	s = currentfe_.invTrJPhiDerMaster(i, 0, iq)*currentfe_.invTrJPhiDerMaster(j, 0, iq);

	   	return s;
	}
};

template <class Integrator, UInt ORDER>
class StiffXY{
  private:

  public:
	inline Real operator() (FiniteElement<Integrator, ORDER>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
	{
	   	Real s = 0;
	   	s = 2*currentfe_.invTrJPhiDerMaster(i, 0, iq)*currentfe_.invTrJPhiDerMaster(j, 1, iq);
	   	return s;
	}
};

template <class Integrator, UInt ORDER>
class StiffYY{
  private:

  public:
	inline Real operator() (FiniteElement<Integrator, ORDER>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
	{
	   	Real s = 0;
	   	s = currentfe_.invTrJPhiDerMaster(i, 1, iq)*currentfe_.invTrJPhiDerMaster(j, 1, iq);

	   	return s;
	}
};




//! A Mass class: a class for the mass operator.
template <class Integrator ,UInt ORDER>
class Mass{
	private:
    //! A reference to FiniteElement<Integrator>
    /*!
     * Stores a reference to the finite element where the mass operator is evaluated.
     */
	//FiniteElement<Integrator, ORDER>& currentfe_;
	public:
	//! A constructor.
	/*!
	 \param fe is a reference to FiniteElement<Integrator>
	 */
	//Mass(FiniteElement<Integrator, ORDER> & fe ):currentfe_(fe){};
    //! A definition of operator () taking two arguments.
    /*!
     * Evaluates the mass operator (i,j) of the current finite elemente.
     * \param i is an unsigned int, current finite element local index
     * \param j is an unsigned int, current finite element local index
     * returns a double.
     */
    inline Real operator() (FiniteElement<Integrator, ORDER>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
    {
    	return currentfe_.phiMaster(i,iq)*  currentfe_.phiMaster(j,iq);
    }
};

//! A vGrad class: a class for the the vectorial Gradient operator.

template<class Integrator, UInt ORDER>
class Grad{
	private:
    //! A reference to FiniteElement<Integrator>
    /*!
     * Stores a reference to the finite element where the stiffness operator is evaluated.
     */
	  //FiniteElement<Integrator, ORDER>& currentfe_;
	public:
	//! A constructor.
	/*!
	 \param fe is a reference to FiniteElement<Integrator>
	 */
	//Grad(FiniteElement<Integrator, ORDER> & fe ):currentfe_(fe){};
    //! A definition of operator () taking three arguments.
    /*!
     * Evaluates the component ic of the vGrad operator (i,j) on the current finite elemente.
     * \param i is an unsigned int, current finite element local index
     * \param j is an unsigned int, current finite element local index
     * \param ic is an unsigned int, vGrad component to be evaluated
     * returns a double.
     */
     inline Real operator() (FiniteElement<Integrator, ORDER>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
     {
    	 return currentfe_.phiMaster(i,iq)*currentfe_.invTrJPhiDerMaster(j,ic,iq);
     }
};


//generic template class wrapper
//! A ETWrapper class: Expression Template Wrapper.
/*!
 * Class that mimic the behaviour of a generic operator defined above: following
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */



template<typename P, typename A>
class EOExpr{
	private:
	  //! "A" is a generic type
	  A a_;
	public:
	//! A constructor.
	/*!
	 * \param object is a constant reference to a generic operator.
	 */
	  EOExpr(const A& a):a_(a){};
	 //! A definition of operator () which takes two arguments.
     /*!
     * Masks the behaviour of the correspondent operator in the above classes.
     * \param i is an unsigned int
     * \param j is an unsigned int
     * returns a P variable.
     */
	  //P operator() (UInt i, UInt j) {return a_(i,j);}
	 //! A definition of operator () which takes three arguments.
     /*!
     * Masks the behaviour of the correspondent operator in the above classes.
     * \param i is an unsigned int
     * \param j is an unsigned int
     * \param ic is an unsigned int
     * returns a P variable.
     */

	  template<typename Integrator, UInt ORDER>
      P operator() (FiniteElement<Integrator, ORDER>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
      {
          return a_(currentfe_, i,j,iq,ic);
      }
	 //! A definition of operator () which takes four arguments.
     /*!
     * Masks the behaviour of the correspondent operator in the above classes.
     * \param i is an unsigned int
     * \param j is an unsigned int
     * \param ic is an unsigned int
     * returns a P variable.
     */
	  //P operator() (UInt i, UInt j, UInt ic1, UInt ic2) { return a_(i,j,ic1,ic2);}
	};

//composition of two wrappers (operator)
//! A ETWBinOp class: Expression Template Wrapper Binary Operation
/*!
 * Class that implements an abstract binary operation defined by Op between two ETWrappers, following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */
template<typename P, typename A, typename B, typename Op>
class EOBinOp{
	private:
		//! "A" is a generic type.
		/*!
		 * Stores the first operand.
		 */
		A a_;
		//! "B" is a generic type.
		/*!
		 * Stores the second operand.
		 */
		B b_;
	public:
	//! A constructor.
	/*!
	 * \param a is a constant reference to a generic type.
	 * \param b is a constant reference to a generic type.
	 */
		EOBinOp(const A& a ,const B& b): a_(a),b_(b){};
	 //! A definition of operator () taking two arguments.
	 /*!
     * \param i is an unsigned int
     * \param j is an unsigned int
     * applies the generic operation defined by the type Op to the two generic objects a_, b_;
     * returns a type P variable
	 */
		template<typename Integrator, UInt ORDER>
	  P operator () (FiniteElement<Integrator, ORDER>& currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
	  {
		  return Op::apply(a_(currentfe_,i,j, iq, ic),b_(currentfe_, i,j, iq, ic));
	  }
	};

template<typename P, class B, class Op>
class EOBinOp<P, Real, B, Op>
{
	Real M_a;
	B M_b;
public:
	EOBinOp(Real a, const B& b):M_a(a),M_b(b) {};

	template<typename Integrator, UInt ORDER>
	inline P operator()(FiniteElement<Integrator, ORDER>& currentfe_, int i, int j, int iq, int ic = 0)
	{
		return Op::apply(M_a,M_b(currentfe_, i, j, iq, ic));
	}
};

template<typename P, class B, class Op>
class EOBinOp<P, Function, B, Op>
{
	const Function& M_a;
	B M_b;
public:
	EOBinOp(const Function& a, const B& b):M_a(a),M_b(b) {};

	template<typename Integrator, UInt ORDER>
	inline P operator()(FiniteElement<Integrator, ORDER>& currentfe_, int i, int j, int iq, int ic = 0)
	{
		UInt globalIndex = currentfe_.getGlobalIndex(iq);
		return Op::apply(M_a(globalIndex),M_b(currentfe_, i, j, iq, ic));
	}
};

//wrappers addition
//! A ETWAdd class: Expression Template Wrapper Addition
/*!
 * Class that defines Addition operation, following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */
template <typename P>
class EOAdd{
	public:
	//! A constructor.
	EOAdd(){}
	//! A stastic inline method taking two arguments.
	/*!
	 *The actual addition operation
	 * \param a is of P type, first addend
	 * \param b is of P type, second addend
	 */
	static inline P apply(P a, P b){ return (a+b); }
};

//multiplication by real scalar
//! A ETWMult class: Expression Template Wrapper Multiplication.
/*!
 * Class that defines Multiplication operation, following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */
template<typename P>
class EOMult{
	public:
	//! A constructor
	EOMult(){}
	 //! A stastic inline method taking two arguments.
	/*!
	 * The actual multiplication operation.
	 * \param a is of P type, first operand
	 * \param b is a Real, second operand
	 */
	  static inline P apply( P a, P b){ return (a*b);}
	  //not needed since in ETRBinOp I did "Op::apply(a_(i,j),b_)"
	  //static inline P apply(Real b, const P&a){ return (a*b);}
	};

//Dot
template<typename P, class A, class B>
class EODotProd
{
	A a_;
	B b_;
	public:
	EODotProd(const A& a, const B& b):a_(a), b_(b) {};
	template<typename Integrator, UInt ORDER>
	inline P operator()(FiniteElement<Integrator, ORDER>& currentfe_, int i, int j, int iq, int ic = 0)
	{
		Real s = 0.;
		for(ic = 0; ic < 2; ++ic)
			s += a_(ic) * b_(currentfe_, i, j, iq, ic);
	return s;
	}
};

//Dot
template<typename P, class B>
class EODotProd<P, Function, B>
{
	const Function& M_a;
	B M_b;
	public:
	EODotProd(const Function& a, const B& b):M_a(a), M_b(b) {};
	template<typename Integrator, UInt ORDER>
	inline Real operator()(FiniteElement<Integrator, ORDER>& currentfe_, int i, int j, int iq, int ic = 0)
	{
		Real s = 0.;
		UInt globalIndex = currentfe_.getGlobalIndex(iq);
		for(ic = 0; ic < 2; ic++)
			s += M_a(globalIndex, ic) * M_b(currentfe_, i, j, iq, ic);
	return s;
	}
};

//operator +
//! Overloading of operator +.
/*!
 * Following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 * Takes two arguments:
 * \param a is const reference ETWrapper<P, A>
 * \param b is const reference ETWrapper<P, A>
 * \return a ETWrapper<P,ETWBinOp<P, ETWrapper<P,A>, ETWrapper<P, B>, ETWAdd<P> > which is resolved at compile time.
 */
template<typename P, typename A, typename B>
EOExpr<P, EOBinOp<P, EOExpr<P, A>, EOExpr<P, B>, EOAdd<P> > >
operator + (const EOExpr<P, A>&  a, const EOExpr<P, B>&  b){

	  typedef EOBinOp<P, EOExpr<P, A>, EOExpr<P, B>, EOAdd<P> > ExprT;
	  return EOExpr<P, ExprT> (ExprT(a,b));
}

template<typename P, typename B>
EOExpr<P, EOBinOp<P, Function, EOExpr<P, B>, EOMult<P> > >
operator * (const Function&  a, const EOExpr<P, B>&  b){

	  typedef EOBinOp<P, Function, EOExpr<P, B>, EOMult<P> > ExprT;
	  return EOExpr<P, ExprT> (ExprT(a,b));
}

template<typename P, typename B>
EOExpr<P, EOBinOp<P, Real, EOExpr<P, B>, EOMult<P> > >
operator * (Real a, const EOExpr<P, B>&  b){

	  typedef EOBinOp<P, Real, EOExpr<P, B>, EOMult<P> > ExprT;
	  return EOExpr<P, ExprT> (ExprT(a,b));
}

template<typename P, typename B>
EOExpr<P, EODotProd<P, Eigen::Matrix<Real,2,1>, EOExpr<P, B> > >
dot(const Eigen::Matrix<Real,2,1>& a, const EOExpr<P, B>&  b){

	  typedef EODotProd<P, Eigen::Matrix<Real,2,1>, EOExpr<P, B> > ExprT;
	  return EOExpr<P, ExprT> (ExprT(a,b));
}

template<typename P, typename B>
EOExpr<P, EODotProd<P, Function, EOExpr<P, B> > >
dot(const Function& a, const EOExpr<P, B>&  b){

	  typedef EODotProd<P, Function, EOExpr<P, B> > ExprT;
	  return EOExpr<P, ExprT> (ExprT(a,b));
}



//!A Assmbler class: discretize a generic differential operator in a sparse matrix
class Assembler{
	private:

	  //MatrixXr oper_mat_;
	  //! A Eigen::sparseMatrix<Real>
	  /*!
	   * Stores the discretized differential operator.
	   */
	  //SpMat    SPoper_mat_;

	public:
	  //! A constructor
	  //Assembler (){};
	  //! A template member taking three arguments: discretize differential operator
	  /*!
	   * \param oper is a template expression : the differential operator to be discretized.
	   * \param mesh is const reference to a MeshHandler<ORDER>: the mesh where we want to discretize the operator.
	   * \param fe is a const reference to a FiniteElement
	   * stores the discretization in SPoper_mat_
	   */

	  //Obsolet
	  //template<UInt ORDER, typename Integrator, typename P, typename A>
	  //void operKernel(EOExpr<P, A> oper,const MeshHandler<ORDER>& mesh,
	  //                   FiniteElement<Integrator, ORDER>& fe );

	  //Return triplets vector
	  template<UInt ORDER, typename Integrator, typename P, typename A>
	  static void operKernel(EOExpr<P, A> oper,const MeshHandler<ORDER>& mesh,
	  	                     FiniteElement<Integrator, ORDER>& fe, SpMat& OpMat);

	  template<UInt ORDER, typename Integrator>
	  static void forcingTerm(const MeshHandler<ORDER>& mesh, FiniteElement<Integrator, ORDER>& fe, const ForcingTerm& u, VectorXr& forcingTerm);

	  //void getOperKernel(UInt i, UInt j){ std::cout<<oper_mat_(i,j)<<std::endl; };
	  //! A normal member taking two arguments: access to element (i,j)
	  /*!
	   *\param i is an unsigned int
       *\param j is an unsigned int
       * print the (i,j) element of SPoper_mat_
	   */
	  //void getSpOperKernel(UInt i, UInt j){ std::cout<<SPoper_mat_.coeffRef(i,j)<<std::endl; };
	  //! A normal member returning a reference to the discretized differential operator
	  /*!
	   * \return a reference to SP_opermat_
	   */
	  //SpMat& getSpOperKernel(){ return (SPoper_mat_); }
	};




#include "matrix_assembler_imp.hpp"

#endif
