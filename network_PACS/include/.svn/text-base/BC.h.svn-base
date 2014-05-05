/*
 * FractureBCHandler.h
 *
 *  Created on: Apr 11, 2011
 *      Author: fumagalli
 */

#ifndef BC_H_
#define BC_H_ 1

#include "Core.h"
#include "UsefulFunctions.h"

class BC
{
public:

    enum
    {
        DIRICHLET_BOUNDARY_NUM = 40, NEUMANN_BOUNDARY_NUM = 50
    };

    BC ( getfem::mesh& mesh,
         const std::string& MeshType,
         const ElementDimension& dimension = MEDIUM );

    inline const sizeVector_Type& getDirichlet ( ) const
    {
        return M_dirichlet;
    }

    inline const size_type& getDirichlet ( const size_type& dof ) const
    {
        return M_dirichlet [ dof ];
    }

    inline const sizeVector_Type& getNeumann ( ) const
    {
        return M_neumann;
    }

    inline const size_type& getNeumann ( const size_type& dof ) const
    {
        return M_neumann [ dof ];
    }

    inline const getfem::mesh_fem& getMeshFEM ( ) const
    {
        return M_meshFEM;
    }

private:

    getfem::mesh_fem M_meshFEM;

    // flags for BC
    sizeVector_Type M_dirichlet;
    // flags for BC
    sizeVector_Type M_neumann;

    // flags for BC
    sizeVector_Type M_extBoundary;
};

typedef BC BC_Type;
typedef boost::shared_ptr<BC> BCPtr_Type;
typedef std::vector<BCPtr_Type> BCPtrContainer_Type;

#endif /* BC_H_ */
