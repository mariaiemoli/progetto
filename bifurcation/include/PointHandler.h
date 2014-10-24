#ifndef _POINTHANDLER_	
#define _POINTHANDLER_

#include <iosfwd>
#include <Eigen/Dense>
#include "Core.h"


typedef Eigen::Vector2d Vector;

const int ndim=2;
 
class PointHandler {
	public:
		static int const myDim=0;
		
		explicit  PointHandler(double x=0, double y=0);// builds a PointHandler
		
		PointHandler(const PointHandler &);
		
		PointHandler & operator=(const PointHandler&);
		
		void setCoordinates(double x, double y); // Sets point coordinates
		
		void getCoordinates(double & x, double & y) const; //Get points coordinates
		
		// We can get the coordinates also by operator []
		double const operator[](int i)const {return M_coor[i];}
		
		double & operator[](int i) {return M_coor[i];}
		
		double const x()const {return M_coor[0];}
		
		double const y()const {return M_coor[1];}
		
		PointHandler operator +=(const PointHandler &);
		
		PointHandler operator -=(const PointHandler &);
		
		friend PointHandler operator +(const PointHandler&, const PointHandler &);
		
		friend PointHandler operator -(const PointHandler&, const PointHandler &);
		
		friend std::ostream & operator <<(std::ostream &, PointHandler const &);
		
		PointHandler operator *(const double &)const;
		
		//! Dot product
		
		double dot(PointHandler const &) const ;
		
		friend PointHandler operator*(const double &, const PointHandler &);
		
		//! Return point as an eigen vector
		Vector asVector()const{return Vector(M_coor[0],M_coor[1]);}

	private:
		double M_coor[ndim];
};



#endif