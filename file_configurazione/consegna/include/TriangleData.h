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

#ifndef __TRIANGLEDATA_H__
#define __TRIANGLEDATA_H__

#include <eigen3/Eigen/Dense>

#include <assert.h>
#include "PointData.h"
#include "StringUtility.h"

/**************************************************************************/
/*  TriangleData.h														  */
/*  Classe che contiene i dati principali per costrire il triangolo       */
/**************************************************************************/

class TriangleData 
{
public:
	static int const myDim=2;
	static const int numVertices=3;
	static const int numSides=3;
		
	// Costruttore nullo
	TriangleData(); 
	
	
	// Costruisce il triangolo a partire da tre punti che rappresentano i vertici
	TriangleData(PointData&,PointData&,PointData&); 
	
	// Costruttore di copia
	TriangleData(const TriangleData&); 
	
	
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
	 * Funzione che calcola l'area del triangolo
	 */
	scalar_type measure() const; 
	
	size_type size () const;
	
	
	/**
	 * Funzione che restituisce uno dei due punti che costituiscono un lato del triangolo
	 * \param size_type edgenum: lato di unteresse
	 * \param size_type endnum: indice del punto che costituisce un vertice 
	 */
	PointData& edgePoint(size_type edgenum, size_type endnum); 
	
	
	/**
	 * Funzione che restituisce uno dei due punti che costituiscono un lato del triangolo
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
	TriangleData& operator=(const TriangleData&);
	
	friend std::ostream & operator <<(std::ostream &, TriangleData const &);
	
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
		
}; // TriangleData


typedef TriangleData TriangleData_Type;												/*!< Classe TriangleData */
typedef std::vector < TriangleData_Type > TriangleDataContainer_Type;				/*!< Vettore di classi TriangleData */
typedef boost::shared_ptr < TriangleData_Type > TriangleDataPtr_Type;				/*!< Puntatore alla classe TriangleData */
typedef std::vector < TriangleDataPtr_Type> TriangleDataPtrContainer_Type;			/*!< Vettore di puntatori alla classe TriangleData */


#endif