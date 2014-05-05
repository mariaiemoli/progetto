/*
 * BCHandler.cc
 *
 *  Created on: Apr 5, 2011
 *      Author: fumagalli
 */

#include "../include/BCHandler.h"

BCHandler::BCHandler ( const BCPtr_Type& mediumBC,
                       const BCPtrContainer_Type& fractureBC ) :
    M_mediumBC(mediumBC), M_fractureBC(fractureBC), M_extBoundaryCut(
            M_fractureBC.size()), M_dirichletCut(M_fractureBC.size()),
            M_neumannCut(M_fractureBC.size())
{
}

// Fix the boundary dividing in cut and uncut
void BCHandler::createBDRegions ( getfem::mesh& mesh )
{
    M_neumannCut.clear();
    M_dirichletCut.clear();

    // Flags for the uncutted region
    M_neumannUncut.clear();
    M_neumannUncut.push_back(NEUMANN_BOUNDARY_UNCUT);
    M_dirichletUncut.clear();
    M_dirichletUncut.push_back(DIRICHLET_BOUNDARY_UNCUT);

    getfem::mesh_region& meshRegionNeumann = mesh.region(
            BC::NEUMANN_BOUNDARY_NUM);
    dal::bit_vector mediumMeshRegionIndex = meshRegionNeumann.index();

    size_type i_cv = 0;
    for ( i_cv << mediumMeshRegionIndex; i_cv != size_type(-1); i_cv
            << mediumMeshRegionIndex )
    {
        for ( size_type jj = 0; jj < 3; ++jj )
        {
            if ( meshRegionNeumann.is_in(i_cv, jj) )
            {
                mesh.region(M_neumannUncut [ 0 ]).add(i_cv, jj);
            }
        }
    }

    getfem::mesh_region& meshRegionDirichlet = mesh.region(
            BC::DIRICHLET_BOUNDARY_NUM);
    mediumMeshRegionIndex = meshRegionDirichlet.index();

    i_cv = 0;
    for ( i_cv << mediumMeshRegionIndex; i_cv != size_type(-1); i_cv
            << mediumMeshRegionIndex )
    {
        for ( size_type jj = 0; jj < 3; ++jj )
        {
            if ( meshRegionDirichlet.is_in(i_cv, jj) )
            {
                mesh.region(M_dirichletUncut [ 0 ]).add(i_cv, jj);
            }
        }
    }
}

void BCHandler::createBDRegionsFractures ( getfem::mesh& mesh )
{
    const size_type numberFractures = M_fractureBC.size();

    getfem::mesh_region& meshRegionNeumann = mesh.region(BC::NEUMANN_BOUNDARY_NUM);
    size_type i_cv;
    dal::bit_vector mediumMeshRegionIndex;

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        mediumMeshRegionIndex = meshRegionNeumann.index();
        i_cv = 0;
        for ( i_cv << mediumMeshRegionIndex; i_cv != size_type(-1); i_cv
                << mediumMeshRegionIndex )
        {
            // If the current element is in the cut region update its boundary conditions
            if ( mesh.region(f + FractureData::FRACTURE).is_in(i_cv) )
            {

                if ( M_neumannCut [ f ].size() <= 0 )
                {
                    M_neumannCut [ f ].push_back(NEUMANN_BOUNDARY_CUT + f);
                }
                for ( size_type jj = 0; jj < 3; ++jj )
                {
                    if ( meshRegionNeumann.is_in(i_cv, jj) )
                    {
                        mesh.region(M_neumannCut [ f ] [ 0 ]).add(i_cv, jj);
                        mesh.region(M_neumannUncut [ 0 ]).sup(i_cv, jj);
                    }
                }
            }
        }
    }

    getfem::mesh_region& meshRegionDirichlet = mesh.region(BC::DIRICHLET_BOUNDARY_NUM);

    for ( size_type f = 0; f < numberFractures; ++f )
    {
        mediumMeshRegionIndex = meshRegionDirichlet.index();
        i_cv = 0;
        for ( i_cv << mediumMeshRegionIndex; i_cv != size_type(-1); i_cv
                << mediumMeshRegionIndex )
        {
            if ( mesh.region(f + FractureData::FRACTURE).is_in(i_cv) )
            {
                if ( M_dirichletCut [ f ].size() <= 0 )
                {
                    M_dirichletCut [ f ].push_back(DIRICHLET_BOUNDARY_CUT + f);
                }
                for ( size_type jj = 0; jj < 3; ++jj )
                {
                    if ( meshRegionDirichlet.is_in(i_cv, jj) )
                    {
                        mesh.region(M_dirichletCut [ f ] [ 0 ]).add(i_cv, jj);
                        mesh.region(M_dirichletUncut [ 0 ]).sup(i_cv, jj);
                    }
                }
            }
        }

    }
}

