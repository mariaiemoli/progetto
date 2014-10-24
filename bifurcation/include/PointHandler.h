#ifndef _POINTHANDLER_H_	
#define _POINTHANDLER_H_

#include <iosfwd>
#include <iostream>
#include <Eigen/Dense>
#include "Core.h"

const size_type ndim=2;
 
class PointHandler {
	public:
		static int const myDim=0;
		
		explicit  PointHandler(scalar_type x=0, scalar_type y=0);// builds a PointHandler
		
		PointHandler(const PointHandler &);
		
		PointHandler & operator=(const PointHandler&);
		
		void setCoordinates(scalar_type x, scalar_type y); // Sets point coordinates
		
		void getCoordinates(scalar_type & x, scalar_type & y) const; //Get points coordinates
		
		// We can get the coordinates also by operator []
		scalar_type const operator[](int i)const {return M_coor[i];}
		
		scalar_type & operator[](int i) {return M_coor[i];}
		
		scalar_type const x()const {return M_coor[0];}
		
		scalar_type const y()const {return M_coor[1];}
		
		PointHandler operator +=(const PointHandler &);
		
		PointHandler operator -=(const PointHandler &);
		
		friend PointHandler operator +(const PointHandler&, const PointHandler &);
		
		friend PointHandler operator -(const PointHandler&, const PointHandler &);
		
		friend std::ostream & operator <<(std::ostream &, PointHandler const &);
		
		PointHandler operator *(const scalar_type &)const;
		
		//! Dot product
		
		scalar_type dot(PointHandler const &) const ;
		
		friend PointHandler operator*(const scalar_type &, const PointHandler &);
		
		//! Return point as an eigen vector
		Vector2d asVector()const{return Vector2d(M_coor[0],M_coor[1]);}

	private:
		scalar_type M_coor[ndim];
};



#endif