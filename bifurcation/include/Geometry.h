/*
 * PROGETTO DI PACS 2014
 *
 * \author Bonomi Claudia
 * 
 * \author Iemoli Maria
 *
 * Problema di Darcy per un network di fratture
 *
 */


#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include "TriangleData.h"
#include "QuadrilaterData.h"
#include "FractureHandler.h"

/**************************************************************************/
/*  Geometry.h															  */
/*  Libreria che costruisce il triangolo 2d dell'intersezione             */
/**************************************************************************/


/**
 * Classe che rappresenta la fine di una frattura
 */ 
class FractureEnd
{
public:
	
	//Costruttore
	FractureEnd( PointData& , scalar_type ); 
	
	inline PointData getPoint ( ) const
	{
		return M_endPoint;
	}
	
	inline scalar_type getThickness ( ) const
	{
		return M_thickness;
	}
	
	//Operatori
	FractureEnd & operator=(const FractureEnd& t)
	{
		if (this !=&t)
		{
			M_endPoint = t.getPoint();
			M_thickness = t.getThickness();
		}

		return *this;
	}

private:
  
	// Punto finale della frattura
	PointData  M_endPoint;
	
	// Spessore della frattura
	scalar_type M_thickness;
		
};

typedef FractureEnd FractureEnd_Type;											/*!< Classe FractureEnd */
typedef std::vector < FractureEnd_Type > FractureEndContainer_Type;				/*!< Vettore di classi FractureEnd */
typedef boost::shared_ptr < FractureEnd_Type > FractureEndPtr_Type;				/*!< Puntatore alla classe FractureEnd */
typedef std::vector < FractureEndPtr_Type> FractureEndPtrContainer_Type;		/*!< Vettore di puntatori alla classe FractureEnd */

/**
 * Classe che rappresenta l'intersezione
 */
class Intersection
{
public:
	
	
	//Costruttori
	/**
	* Le fratture devono essere date in senso orario
	*/ 
	Intersection( FractureEndContainer_Type const & fractureEnd, PointData const & intersectionPoint);
	
	Intersection()
	{};
	
	/**
	 * Funzione che costruisce il triangolo di intersezione partendo dal vettore delle fratture che si intersecano
	 */
	void setIntersection( FracturePtrContainer_Type& M_FracturesSet );
	
	/**
	 * Funzione che costruisce il triangolo di intersezione partendo dal vettore delle fratture che si intersecano
	 */
	void setTriangleIntersection( FracturePtrContainer_Type& M_FracturesSet );


	/**
	 * Funzione che costruisce il quadrilatero di intersezione partendo dal vettore delle fratture che si intersecano
	 */
	void setQuadrilaterIntersection( FracturePtrContainer_Type& M_FracturesSet );


	/**
	 * Funzione che costruisce il triangolo di intersezione 
	 */
	TriangleData const & computeIntersectionTriangle();
	

	/**
	 * Funzione che costruisce il quadrilatero di intersezione 
	 */
	QuadrilaterData const & computeIntersectionQuadrilater();
	
	
	/**
	 * Funzione che costruisce il triangolo di intersezione a partire da un triangolo già esistente
	 */
	void setTriangle ( const TriangleData& triangle)
	{
		M_intersectionTriangle = triangle;
		
		return;
	}
	
	
	/**
	 * Funzione che restituisce il triangolo di intersezione
	 */
	TriangleData const & intersectionTriangle()const 
	{
		return M_intersectionTriangle;
	}
	
	/**
	 * Funzione che restituisce il triangolo di intersezione
	 */
	QuadrilaterData const & intersectionQuadrilater()const 
	{
		return M_intersectionQuadrilater;
	}
	
	/**
	 * Funzione che restituisce il punto di intersezione
	 */
	PointData getPointIntersection( ) const
	{
		return M_intersection;
	}
	
	
	Vector2d getNormal( const size_type i ) const
	{
		return M_normals [ i ];
	}
	
	Vector2d getTangent( const size_type i ) const
	{
		return M_tangents [ i ];
	}
	
	/**
	 * Funzione che ordina le fratture in senso orario
	 */
	void sortFractures ( FracturePtrContainer_Type& M_FracturesSet );
	
	
	/**
	 * Funzione che mi restituisce vero o falso a seconda se una frattura si trovi o meno nella parte positiva di un'altra
	 */
	bool isPos ( FractureHandlerPtr_Type& fracture, FractureHandlerPtr_Type& otherfracture);
	
	Intersection & operator =(const Intersection & Int)
	{
		this-> M_fractures = Int.M_fractures;
		this-> M_intersection = Int.M_intersection;
		
		this-> M_tangents [ 0 ] = Int.M_tangents [ 0 ];
		this-> M_tangents [ 1 ] = Int.M_tangents [ 1 ];
		this-> M_tangents [ 2 ] = Int.M_tangents [ 2 ];
		this-> M_tangents [ 3 ] = Int.M_tangents [ 3 ];
		
		this-> M_normals [ 0 ] = Int.M_normals [ 0 ];
		this-> M_normals [ 1 ] = Int.M_normals [ 1 ];
		this-> M_normals [ 2 ] = Int.M_normals [ 2 ];
		this-> M_normals [ 3 ] = Int.M_normals [ 3 ];
		
		this-> M_intersectionTriangle = Int.M_intersectionTriangle;
		this-> M_intersectionQuadrilater = Int.M_intersectionQuadrilater;
		
		return *this;
		
	}
  

private:
	  
	// Insieme delle fratture che si interscano
	FractureEndContainer_Type M_fractures;
	
	// Punto di intersezione
	PointData M_intersection;
	
	Vector2d M_tangents [ 4 ];
	Vector2d M_normals [ 4 ];
	
	// Triangolo di intersezione
	TriangleData M_intersectionTriangle;
	
	// Quadrilatero di intersezione
	QuadrilaterData M_intersectionQuadrilater;

};

typedef Intersection Intersection_Type;											/*!< Classe Intersection */
typedef std::vector < Intersection_Type > IntersectionContainer_Type;			/*!< Vettore di classi Intersection */
typedef boost::shared_ptr < Intersection_Type > IntersectionPtr_Type;			/*!< Puntatore alla classe Intersection */
typedef std::vector < IntersectionPtr_Type> IntersectionPtrContainer_Type;		/*!< Vettore di puntatori alla classe Intersection */

#endif		/* __GEOMETRY_H__ */