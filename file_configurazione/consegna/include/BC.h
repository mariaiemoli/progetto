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

#ifndef BC_H_
#define BC_H_ 1

#include "Core.h"

/**************************************************************************/
/*  BC.h															      */
/*  Classe che introduce le condizioni al bordo sul problema              */
/**************************************************************************/

class BC
{
public:

	//static sizeVector_Type DEFAULT_VECTOR;

    enum
    {
        DIRICHLET_BOUNDARY_NUM = 40,
        NEUMANN_BOUNDARY_NUM = 50
    };


    // Costruttore
    BC ( getfem::mesh& mesh,
         const std::string& MeshType,
		 const sizeVector_Type DOFs /*= DEFAULT_VECTOR*/ );


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
    sizeVector_Type M_neumann;
    sizeVector_Type M_extBoundary;
};

typedef BC BC_Type;												/*!< Classe BC */
typedef boost::shared_ptr<BC> BCPtr_Type;						/*!< Puntatore alla classe BC */
typedef std::vector<BCPtr_Type> BCPtrContainer_Type;			/*!< Vettore di puntatori alla classe BC */

#endif /* BC_H_ */
