#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>

class Point
{

public:
   Point ( const double& x, const double& y):
   M_x ( x ),
   M_y ( y )
   {}

   inline double x() const { return M_x; }
   inline double y() const { return M_y; }

private:
    double M_x;
    double M_y;
};

class Triangle
{
public:
    Triangle ( const unsigned int& id1, const unsigned int& id2, const unsigned int& id3 )
    {   M_id.resize ( 3 );
        M_id[0] = id1; M_id[1] = id2; M_id[2] = id3; }

    inline unsigned int id ( const unsigned int& i ) const { return M_id.at( i ); }

private:
    std::vector<unsigned int> M_id; 
};

typedef std::vector<Point> Points;
typedef std::vector<Triangle> Triangles;

void readMeshFreeFem( const std::string& fileMesh, Points& points, Triangles& triangles );

void writeMeshGetFem( const std::string& fileMesh, const Points& points, const Triangles& triangles );

int main( int argc, char* argv[] )
{

    if ( argc < 3 )
    {
        std::cerr << "More arguments are needed" << std::endl;
        exit(1);
    }
    
    const std::string fileFreeFem = argv[1];
    const std::string fileMesh = argv[2];

    std::cout << "Call FreeFem++ for generate the mesh" << std::endl;
    system( ( "FreeFem++ " + fileFreeFem ).data () );

    Points points;
    Triangles triangles;
  
    std::cout << "Read data from FreeFem++ mesh file" << std::endl;
    readMeshFreeFem ( fileMesh, points, triangles );

    std::cout << "Write data to GetFem mesh file" << std::endl;
    writeMeshGetFem ( fileMesh, points, triangles );

    return 0;
} // main

void readMeshFreeFem ( const std::string& fileMesh, Points& points, Triangles& triangles )
{
    std::ifstream file;
    file.open ( (fileMesh + ".msh" ).data (), std::ifstream::in );

    unsigned int nbPoints;
    file >> nbPoints;
    points.reserve ( nbPoints );

    unsigned int nbTriangles;
    file >> nbTriangles;
    triangles.reserve ( nbTriangles );

    double dummy;
    file >> dummy; 

    std::cout << "Number of points " << nbPoints << std::endl << " Number of triangles " << nbTriangles << std::endl;

    double x, y;
    for ( unsigned int i = 0; i < nbPoints; ++i )
    {
        file >> x >> y >> dummy;
        points.push_back ( Point ( x, y ) );
    }

    unsigned int id1, id2, id3;
    for ( unsigned int i = 0; i < nbTriangles; ++i )
    {
        file >> id1 >> id2 >> id3 >> dummy;
        triangles.push_back ( Triangle ( id1, id2, id3 ) );
    }

    file.close ();
} // readMeshFreeFem

void writeMeshGetFem ( const std::string& fileMesh, const Points& points, const Triangles& triangles )
{
    std::ofstream file;
    file.open ( ( fileMesh ).data () );

    file << "% GETFEM MESH FILE" << std::endl << "% GETFEM VERSION 4.1.1" << std::endl;
    file << std::endl << std::endl << std::endl;

    file << "BEGIN POINTS LIST" << std::endl << std::endl;

    for ( unsigned int i = 0; i < points.size(); ++i ) 
    {
        file << "  POINT  " << i << " " << points[i].x() << " " << points[i].y() << std::endl;
    }

    file << std::endl << "END POINTS LIST" << std::endl;
    file  << std::endl << std::endl << std::endl;

    file << "BEGIN MESH STRUCTURE DESCRIPTION" << std::endl << std::endl;

    for ( unsigned int i = 0; i < triangles.size(); ++i )
    {
        file << "CONVEX " << i << " 'GT_PK(2,1)' " << triangles[i].id(0) - 1 << " " 
                                                   << triangles[i].id(1) - 1 << " "
                                                   << triangles[i].id(2) - 1 << std::endl;
    }

    file << std::endl << "END MESH STRUCTURE DESCRIPTION";

    file.close ();
} // writeMeshGetFem
