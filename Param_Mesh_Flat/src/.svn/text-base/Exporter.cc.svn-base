/*
 * Exporter.cc
 *
 *  Created on: Apr 13, 2011
 *      Author: fumagalli
 */

#include "../include/Exporter.h"

Exporter::
Exporter ( const GetPot& dataFile, const std::string& section ) :
    M_vtkFolder(dataFile((section + "folderVTK").data(), "./vtk/"))
{
}

void Exporter::
spy ( const sparseMatrixPtr_Type& matrix, const std::string& nameFile ) const
{
    gmm::MatrixMarket_IO::write(nameFile.c_str(), *matrix);
} // spy

void Exporter::
spy ( const scalarVectorPtr_Type& vector, const std::string& nameFile ) const
{
    std::ofstream file;
    file.open(nameFile.c_str());
    for ( size_type i = 0; i < vector->size(); ++i )
    {
        file << std::setprecision(15) << vector->at(i) << std::endl;
    }
    file.close();
} // spy

void Exporter::
meshRegion ( const getfem::mesh& mesh, const std::string& nameFile ) const
{
        // set up a P0 space, each elements contains the region mesh value
        GFMeshFEM_Type meshFEM ( mesh );
        getfem::pfem FEType = getfem::fem_descriptor ( "FEM_PK(2,0)" );
        meshFEM.set_finite_element ( FEType );

        // vector to store the values of the region mesh
        const size_type nbElements = meshFEM.nb_basic_dof();
        scalarVector_Type regionMesh ( nbElements, 0 );

        // Extract the regions flags from the mesh
        dal::bit_vector regionsIndexBitVector ( mesh.regions_index() );
        sizeVector_Type regionsIndex;
        fromBitVectorToStdVector ( regionsIndexBitVector, regionsIndex );

        for ( size_type i = 0; i < regionsIndex.size(); ++i )
        {
                getfem::mesh_region meshCurrentRegion ( mesh.region ( regionsIndex [i] ) );
                // check if the current region is made just with convexes, avoid boundary regions
                if ( meshCurrentRegion.is_only_convexes() == true )
                {
                        dal::bit_vector currentRegionBitVector ( meshCurrentRegion.index() );
                        size_type i_cv = 0;
                        for ( i_cv << currentRegionBitVector; i_cv != bgeot::size_type(-1); i_cv << currentRegionBitVector )
                        {
                                const size_type position = (meshFEM.ind_basic_dof_of_element( i_cv ))[0];
                                regionMesh [ position ] += regionsIndex [i];
                        }
                }
        }

        // Export as a solution
        exportSolution ( M_vtkFolder + nameFile, "RegionMesh", meshFEM, regionMesh );

} // meshRegion
