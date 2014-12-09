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


#ifndef LEVELSETHANDLER_H_
#define LEVELSETHANDLER_H_ 1

#include "LevelSetData.h"


/**************************************************************************/
/*  LevelSetHandler.h													  */
/*  Classe che inizializza e gestisce il levelset che rappresenta una  	  */
/*  frattura                 			  								  */
/**************************************************************************/

class LevelSetHandler
{
public:

    LevelSetHandler ( const GetPot& dataFile, const std::string& section =  "fractureData/", 
    				  const std::string& sectionLevelSet = "levelSet/" );

    
    /**
     * Funzione che inizializza il level set. Partendo dalle informazioni contenute nel file data costrusce le mesh e i metodi di integrazione
     * \param getfem::mesh& mediumMesh: mesh di supporto del mezzo
     * \param std::string& mediumIntegrationTypeVelocity: metodo di integrazione per la velocità 
     * \param getfem::mesh_fem& mediumMeshFEMPressure: mesh di integrazione per la pressione
     * \param getfem::mesh_fem& mediumMeshFEMVelocity: mesh di integrazione per la velocità
     */
    void init ( getfem::mesh& mediumMesh,
			    const std::string& mediumIntegrationTypeVelocity,
			    const getfem::mesh_fem& mediumMeshFEMPressure,
			    const getfem::mesh_fem& mediumMeshFEMVelocity );


    inline const GFMeshLevelSet_Type& getMesh ( ) const
    {
        return *M_mesh;
    }


    inline GFMeshLevelSet_Type& getMesh ( )
    {
        return *M_mesh;
    }


    inline const GFIntegrationMethodLevelSet_Type& getIntegrationMethodInside ( ) const
    {
        return *M_integrationMethodInside;
    }


    inline const GFIntegrationMethodLevelSet_Type& getIntegrationMethodOutside ( ) const
    {
        return *M_integrationMethodOutside;
    }


    inline const GFLevelSet_Type& getLevelSet ( ) const
    {
        return *M_levelSet;
    }


    inline GFLevelSet_Type& getLevelSet ( )
    {
        return *M_levelSet;
    }


    inline const GFIntegrationMethodLevelSet_Type& getIntegrationMethod ( ) const
    {
        return *M_integrationMethod;
    }


    inline const scalarVector_Type& getBaricenterValue ( ) const
    {
        return M_baricenterValue;
    }


    inline const scalar_type& getBaricenterValue ( const size_type& dof ) const
    {
        return M_baricenterValue [ dof ];
    }


    inline const scalarVector_Type& getDOFValue ( ) const
    {
        return M_DOFValue;
    }


    inline const scalar_type& getDOFValue ( const size_type& dof ) const
    {
        return M_DOFValue [ dof ];
    }


    inline LevelSetDataPtr_Type& getData ( )
    {
        return M_data;
    }


private:

    // Data
    LevelSetDataPtr_Type M_data;

    // M_mediumMesh, level set
    GFMeshLevelSetPtr_Type M_mesh;
    // Il level set che definisce la frattura
    GFLevelSetPtr_Type M_levelSet;

    // valore del level set nei baricentri
    scalarVector_Type M_baricenterValue;
    // valore del level set nei gradi di libertà di v
    scalarVector_Type M_DOFValue;

    // metodo di integrazione per il level set sul bordo
    GFIntegrationMethodLevelSetPtr_Type M_integrationMethod;
    // metodo di integrazione per il level set " interno "
    GFIntegrationMethodLevelSetPtr_Type M_integrationMethodInside;
    // metodo di integrazione per il level set " esterno "
    GFIntegrationMethodLevelSetPtr_Type M_integrationMethodOutside;

};

typedef LevelSetHandler LevelSetHandler_Type;									/*!< Classe LevelSetHandler */
typedef boost::shared_ptr<LevelSetHandler_Type> LevelSetHandlerPtr_Type;		/*!< Puntatore alla classe LeverSetHandler */

#endif /* LEVELSETHANDLER_H_ */
