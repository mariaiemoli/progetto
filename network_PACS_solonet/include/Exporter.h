/*
 * Exporter.h
 *
 *  Created on: Apr 13, 2011
 *      Author: fumagalli
 */

#ifndef EXPORTER_H_
#define EXPORTER_H_ 1

#include "Core.h"
#include "MeshHandler.h"
#include <fstream>

class Exporter
{
public:

	// costruttore
    Exporter ( const GetPot& dataFile, const std::string& section = "" );

    inline const std::string& getFolder ( ) const
    {
        return M_vtkFolder;
    }

    void	// esporto la matrice
            spy ( const sparseMatrixPtr_Type& matrix,
                  const std::string& nameFile ) const;

    void	// esporto il vettore
            spy ( const scalarVectorPtr_Type& vector,
                  const std::string& nameFile ) const;

    void 	// esporto il "pezzo" di mesh che mi interessa
    		meshRegion ( const getfem::mesh& mesh,
                      const std::string& nameFile = "RegionMesh.vtk" ) const;

private:

    const std::string M_vtkFolder;

};

typedef Exporter Exporter_Type;
typedef boost::shared_ptr<Exporter_Type> ExporterPtr_Type;

#endif /* EXPORTER_H_ */
