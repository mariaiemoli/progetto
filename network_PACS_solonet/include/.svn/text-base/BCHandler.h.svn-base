/*
 * BCHandler.h
 *
 *  Created on: Apr 5, 2011
 *      Author: fumagalli
 */

#ifndef BCHANDLER_H_
#define BCHANDLER_H_ 1

#include "Core.h"
#include "BC.h"
#include "FractureData.h"

class BCHandler
{
public:

    enum
    {
        DIRICHLET_BOUNDARY_CUT = 400,
        NEUMANN_BOUNDARY_CUT = 500,
        DIRICHLET_BOUNDARY_UNCUT = 4000,
        NEUMANN_BOUNDARY_UNCUT = 5000
    };

    BCHandler ( const BCPtr_Type& mediumBC,
                const BCPtrContainer_Type& fractureBC );

    void createBDRegions ( getfem::mesh& mesh );

    void createBDRegionsFractures ( getfem::mesh& mesh );

    inline const BCPtr_Type& getMediumBC ( ) const
    {
        return M_mediumBC;
    }

    inline const BCPtr_Type& getFractureBC ( const size_type& id ) const
    {
        return M_fractureBC [ id ];
    }

    inline const sizeVector_Type& getDirichletUncut ( ) const
    {
        return M_dirichletUncut;
    }

    inline const size_type& getDirichletUncut ( const size_type& dof ) const
    {
        return M_dirichletUncut [ dof ];
    }

    inline const sizeVector_Type& getNeumannUncut ( ) const
    {
        return M_neumannUncut;
    }

    inline const size_type& getNeumannUncut ( const size_type& dof ) const
    {
        return M_neumannUncut [ dof ];
    }

    inline const sizeVector_Type& getDirichletCut ( const size_type& id ) const
    {
        return M_dirichletCut [ id ];
    }

    inline const size_type& getDirichletCut ( const size_type& id,
                                              const size_type& dof ) const
    {
        return M_dirichletCut [ id ] [ dof ];
    }

    inline const sizeVector_Type& getNeumannCut ( const size_type& id ) const
    {
        return M_neumannCut [ id ];
    }

    inline const size_type& getNeumannCut ( const size_type& id,
                                            const size_type& dof ) const
    {
        return M_neumannCut [ id ] [ dof ];
    }

private:

    BCPtr_Type M_mediumBC;
    BCPtrContainer_Type M_fractureBC;

    // flags for BC  bordo tagliato
    sizeVector_Type M_extBoundaryCut;

    // flags for BC bordo non tagliato
    sizeVector_Type M_extBoundaryUncut;

    // flags for BC   divido in parte tagliata e non
    sizeVectorContainer_Type M_dirichletCut;
    // flags for BC
    sizeVectorContainer_Type M_neumannCut;

    // flags for BC
    sizeVector_Type M_dirichletUncut;
    // flags for BC
    sizeVector_Type M_neumannUncut;

};

typedef BCHandler BCHandler_Type;
typedef boost::shared_ptr<BCHandler_Type> BCHandlerPtr_Type;

#endif /* BCHANDLER_H_ */
