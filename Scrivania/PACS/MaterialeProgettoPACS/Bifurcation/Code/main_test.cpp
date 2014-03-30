#include <iostream>
#include "biforc.hpp"


int main ()
{
  using namespace std;
  using namespace Darcy;		// definito in biforc.hpp
  using namespace Geometry;		// definito in geo.hpp
  
  
  // dichiaro i punti medi dei lati del triangolo che delimita la biforcazione
  Point e0(0.,-1.);			
  Point e1(-1.,0.);
  Point e2(1.,0.);
  Point ip(0.,0.);	// baricentro del triangolo
  
  
  // dichiaro dove termina la frattura, la delimito coi punti medi dei lati del triangolo e l'ampiezza di ogni canale
  // FractureEnd è una struct con due campi, la classe Point e un double per lo spessore, in pratica determino i lati del triangolo
  FractureEnd f0{e0,0.1};
  FractureEnd f1{e1,0.1};
  FractureEnd f2{e2,0.05};
  
  
  // classe che definisce l'intersezione a partire dai lati della frattura e dal baricentro del triangolo
  Intersection intersection(f0,f1,f2,ip);
  
  
  // costruisco il triangolo che definisce l'intersezione
  Triangle triangle=intersection.computeIntersectionTriangle();
  std::cout<<triangle<<std::endl;


  // definisco la matrice di permeabilità, matrice 2x2
  Matrix2d permeability;

  //  permeability<<2500.3,2499.8, 
  //            2499.8, 2500.3;
  //  permeability<<2500.0,0.0, 
  //            0.0, 0.2;
  permeability<<1000.0,0.0, 
                0.0, 1.0;
  cout<<"**** Eigenvalues of permeability matrix"<<endl;
  cout<<permeability.eigenvalues()<<endl;


  // punti di intersezione, vertici del triangolo
  Point p0(0.0,0.0);
  Point p1(1.0,0.0);
  Point p2(0,1.0);

  // triangle.setPoint(0,p0);
  //triangle.setPoint(1,p1);
  //triangle.setPoint(2,p2);

  
  Biforc biforc(permeability);
  biforc.setTriangle(triangle);
  biforc.computeT();
  cout<<" ************ Permeability **************"<<endl;
  cout<<biforc.K()<<endl;
  cout<<" ************ Matrix C     **************"<<endl;
  cout<<biforc.C()<<endl;
  cout<<" ************ Matrix N **************"<<endl;
  cout<<biforc.N()<<endl;
  cout<<" ************ Matrix Qc **************"<<endl;
  cout<<biforc.Qc()<<endl;
  cout<<" ************ Matrix Pc **************"<<endl;
  cout<<biforc.Pc()<<endl;
  cout<<" ************ Matrix T  **************"<<endl;
  cout<<biforc.T()<<endl; 
  cout<<" ************eigenvalues  **************"<<endl;
  cout<<biforc.T().eigenvalues()<<endl; 
  cout<<" ************pressure coeff  **************"<<endl;
  cout<<pressureCoeff(biforc.T())<<endl; 

  
}
