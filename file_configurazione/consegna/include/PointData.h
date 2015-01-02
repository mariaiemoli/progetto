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

#ifndef _POINTHANDLER_H_	
#define _POINTHANDLER_H_

#include <iosfwd>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "Core.h"


/**************************************************************************/
/*  PointData.h															  */
/*  Classe che rappresenta un punto geometrico			                  */
/**************************************************************************/
 
class PointData {
	public:
	
		static int const myDim=0;
		
		//Costruttori
		explicit  PointData(scalar_type x=0, scalar_type y=0);
		
		PointData(const PointData &);
		
		//Setta le coordinate del punto
		void setCoordinates(scalar_type x, scalar_type y); 
		
		// Operatori e metodi che mi permettono di avere le coordinate del punto
		void getCoordinates(scalar_type & x, scalar_type & y) const; 
		
		scalar_type const operator[](int i)const 
			{
				return M_coor[i];
			} //operator []
		
		scalar_type & operator[](int i) 
			{
				return M_coor[i];
			} //operator []
		
		
		scalar_type const x()const 
			{
				return M_coor[0];
			} //x()
		
		scalar_type const y()const 
			{
				return M_coor[1];
			} //y()
		
		//Operatori
		PointData & operator=(const PointData&);
		
		PointData operator +=(const PointData &);
		
		PointData operator -=(const PointData &);
		
		friend PointData operator +(const PointData&, const PointData &);
		
		friend PointData operator -(const PointData&, const PointData &);
		
		friend std::ostream & operator <<(std::ostream &, PointData const &);
		
		PointData operator *(const scalar_type &)const;
		
		//! Dot product
		
		scalar_type dot(PointData const &) const ;
		
		friend PointData operator*(const scalar_type &, const PointData &);
		
		//! Conversione Point da vettore di scalari a Vector2d eigen
		Vector2d asVector()const
			{ 
				return Vector2d(M_coor[0],M_coor[1]);
			}//asVector

	private:
		
		//const size_type ndim;
		scalarVector_Type M_coor;
};

typedef PointData PointData_Type;
typedef std::vector < PointData_Type > PointDataContainer_Type;
typedef boost::shared_ptr < PointData_Type > PointDataPtr_Type;
typedef std::vector < PointDataPtr_Type> PointDataPtrContainer_Type;


#endif