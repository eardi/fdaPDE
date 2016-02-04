#ifndef MESH_H_
#define MESH_H_

#include "fdaPDE.h"
#include "mesh_objects.h"


using std::vector;

//!  This class gives an object-oriented reading interface to the output of the library Triangle (Jonathan Richard Shewchuk).
/*!
 * The template parameters specify the order of its elemnts.
 * The aim of this class is to do not introduce any initialization overhead,
 * beacuse it will be called many time during the execution of a R script
*/
template <UInt ORDER>
class MeshHandler {
public:
	typedef int UInt;
	//! A constructor.
    /*!
      * The constructor permits the initialization of the mesh from an R object
      * constructed with the TriLibrary (our R wrapper for the Triangle library)
    */
    
    MeshHandler(Real* points, UInt* edges, UInt* triangles, UInt* neighbors, UInt num_nodes, UInt num_edges, UInt num_triangles):
			points_(points), edges_(edges), triangles_(triangles), neighbors_(neighbors), num_nodes_(num_nodes), num_edges_(num_edges), num_triangles_(num_triangles) {};
    
    #ifdef R_VERSION_
	MeshHandler(SEXP Rmesh);
	#endif
	
	~MeshHandler(){};
	
	//! A normal member returning an unsigned integer value.
    /*!
      \return The number of nodes in the mesh
    */    
    UInt num_nodes() const {return num_nodes_;}

	//! A normal member returning an unsigned integer value.
    /*!
      \return The number of nodes in the mesh
    */    
    UInt num_triangles() const {return num_triangles_;}
    
    //! A normal member returning an unsigned integer value.
    /*!
      \return The number of edges in the mesh
    */ 
    UInt num_edges() const {return num_edges_;}
    
    //! A normal member returning a Point
    /*!
     * \param id an Id argument 
      \return The point with the specified id
    */ 
    Point getPoint(Id id);
    
    //! A normal member returning an Edge
    /*!
     * \param id an Id argument 
      \return The edge with the specified id
    */ 
    Edge getEdge(Id id);
    
    //! A normal member returning a Triangle
    /*!
     * \param id an Id argument 
      \return The triangle with order coerent to that of the mesh with the specified id
    */ 
    Triangle<ORDER * 3>  getTriangle(Id id) const;
    
    //The "number" neighbor of triangle i is opposite the "number" corner of triangle i
    //! A normal member returning the Neighbors of a triangle
    /*!
     * \param id the id of the triangle
     * \param number the number of the vertex
      \return The triangle that has as an edge the one opposite to the specified
      vertex
    */ 
    Triangle<ORDER * 3> getNeighbors(Id id_triangle, UInt number) const;
     
    void printPoints(std::ostream & out);
    void printEdges(std::ostream & out);
    void printTriangles(std::ostream & out);
    void printNeighbors(std::ostream & out);
    
     //! A normal member returning the triangle on which a point is located
    /*!
     * This method implements a simply research between all triangle of the mesh
     * \param point the point we want to locate
      \return The triangle that contains the point
    */ 
    Triangle<ORDER * 3> findLocationNaive(Point point) const;
    
     //! A normal member returning the triangle on which a point is located
    /*!
     * This method implements a Visibility Walk Algorithm (further details in: Walking in a triangulation, Devillers et al)
     * \param point the point we want to locate
     * \param starting_triangles a vector of points that specifies the poposed starting 
     * points for the walking algorithm
      \return The triangle that contains the point
    */ 
    Triangle<ORDER * 3> findLocationWalking(const Point& point, const Triangle<ORDER * 3>& starting_triangle) const;
    
    //int readMesh(std::string const & file);
	//double measure()const;
	//bool checkmesh()const;
private:
	#ifdef R_VERSION_
	SEXP mesh_;
	#endif
	Real *points_;
	UInt *edges_;
	UInt *triangles_;
	UInt *neighbors_;
	
	UInt *border_edges; //contiene lista id_edges al bordo
	UInt num_nodes_, num_edges_, num_triangles_;
	
};

#include "mesh_imp.h"

#endif
