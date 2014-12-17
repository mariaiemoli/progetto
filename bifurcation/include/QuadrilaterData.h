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


#ifndef __QUADRILATERDATA_H__
#define __QUADRILATERDATA_H__ 

#include <eigen3/Eigen/Dense>

#include <assert.h>
#include <math.h>
#include "PointData.h"
#include "StringUtility.h"


/**************************************************************************/
/*  QuadrilaterData.h											  		  */
/*  Classe che contiene i dati principali per costrire il quadrilatero 	  */
/*  di intersezione                 									  */
/**************************************************************************/

class QuadrilaterData
{
public:
	
	static int const myDim=2;
	static const int numVertices=4;
	static const int numSides=4;
		
	// Costruttore nullo
	QuadrilaterData(); 
	
	// Costruisce il quadrilatero a partire da quattro punti che rappresentano i vertici
	QuadrilaterData( PointData&, PointData&, PointData&, PointData& ); 
	
	// Costruttore di copia
	QuadrilaterData( const QuadrilaterData& ); 

	/**
	 * Funzione cbe mi permette di accedere al punto i-esimo del triangolo
	 */
	PointData& getPoint( size_type i ) 
	{
		assert ( i < M_point.size() );
		
		return M_point[i];
	}
		
		
	/**
	 * Funzione che setta un punto del triangolo
	 * \param size_type i: indice del punto da modificare
	 * \param PointData const &  p: punto da inserire
	 */
	void setPoint(size_type i, PointData const &  p);

	/**
	 * Funzione che calcola l'area del quadrilatero
	 */
	scalar_type measure() const; 

	
	/**
	 * Funzione che restituisce uno dei due punti che costituiscono un lato del quadrilatero
	 * \param size_type edgenum: lato di unteresse
	 * \param size_type endnum: indice del punto che costituisce un vertice 
	 */
	PointData& edgePoint(size_type edgenum, size_type endnum); 
	
	
	/**
	 * Funzione che restituisce uno dei due punti che costituiscono un lato del quadrilatero
	 * \param size_type edgenum: lato di unteresse
	 * \param size_type endnum: indice del punto che costituisce un vertice 
	 */
	PointData const & edgePoint(size_type edgenum,size_type endnum) const;
	
		
	/**
	 * Funzione che restituisce le coordinate del baricentro del triangolo
	 */
	PointData baricenter()const;


	// Numerazione dei lati
	static size_type edge(size_type edgeNum, size_type endNum);
	
	
	/**
	 * Funzione che calcola il punto medio di un lato
	 */
	PointData edgeMean(size_type edgeNum) const;
	
	
	/**
	 * Funzione che restituisce il vettore che collega il punto medio di un lato con il baricentro
	 */
	Vector2d c(size_type edgeNum) const;

	
	/**
	 * Funzione che restituisce il vettore che rappresenta la normale uscente da un lato
	 */
	Vector2d unscaledNormal(size_type edgeNum)const;
	
	//Operatori
	QuadrilaterData& operator=(const QuadrilaterData&);
	
	friend std::ostream & operator <<(std::ostream &, QuadrilaterData const &);
	
	
	// Can be used ONLY if empty()==false
	PointData const & operator[](size_type i) const 
	{
		return M_point[i];
	}
	
	// Can be used ONLY if empty()==false
	PointData & operator[](size_type i)
	{
		return M_point[i];
	}

	
private:
	
	// Vettore dei punti che costituiscono il triangolo
	PointDataContainer_Type  M_point;
	
	// Lati del triangolo
	static size_type const M_edge[numSides][2];
	
};



typedef QuadrilaterData QuadrilaterData_Type;												/*!< Classe QuadrilaterData */
typedef std::vector < QuadrilaterData_Type > QuadrilaterDataContainer_Type;					/*!< Vettore di classi QuadrilaterData */
typedef boost::shared_ptr < QuadrilaterData_Type > QuadrilaterDataPtr_Type;					/*!< Puntatore alla classe QuadrilaterData */
typedef std::vector < QuadrilaterDataPtr_Type> QuadrilaterDataPtrContainer_Type;			/*!< Vettore di puntatori alla classe QuadrilaterData */


#endif