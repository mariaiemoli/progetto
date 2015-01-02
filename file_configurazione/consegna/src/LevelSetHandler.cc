
#include "../include/LevelSetHandler.h"

/**************************************************************************/
/*  LevelSetHandler.cc													  */
/*  Classe che inizializza e gestisce il levelset che rappresenta una  	  */
/*  frattura                 			  								  */
/**************************************************************************/

LevelSetHandler::LevelSetHandler ( const GetPot& dataFile,
                                   const std::string& section,
                                   const std::string& sectionLevelSet ) :
                                   M_data( new LevelSetData_Type ( dataFile, section, sectionLevelSet ) )
{}



void LevelSetHandler::init ( getfem::mesh& mediumMesh,
                             const std::string& mediumIntegrationTypeVelocity,
                             const getfem::mesh_fem& mediumMeshFEMPressure,
                             const getfem::mesh_fem& mediumMeshFEMVelocity )
{
    M_mesh.reset ( new GFMeshLevelSet_Type ( mediumMesh ) );

    M_levelSet.reset ( new GFLevelSet_Type ( mediumMesh ) );

    M_integrationMethod.reset ( new GFIntegrationMethodLevelSet_Type ( *M_mesh, getfem::mesh_im_level_set::INTEGRATE_BOUNDARY ) );

    M_integrationMethodInside.reset ( new GFIntegrationMethodLevelSet_Type( *M_mesh, getfem::mesh_im_level_set::INTEGRATE_INSIDE ) );

    M_integrationMethodOutside.reset ( new GFIntegrationMethodLevelSet_Type ( *M_mesh, getfem::mesh_im_level_set::INTEGRATE_OUTSIDE ) );

    // Level sets

    // Descrizione del metodo di integrazione per la funzione level set
    getfem::pintegration_method integrationType = getfem::int_method_descriptor( M_data->getIntegrationTypeSimplex() );

    // Definiamo che la M_mediumMesh M_levelSetMesh è tagliata da un level set
    M_mesh->add_level_set( *M_levelSet );

    // Descrizione del metodo di integrazione per la velocità
    getfem::pintegration_method mediumIntegrationVelocity = getfem::int_method_descriptor( mediumIntegrationTypeVelocity );

    // Descrizione del metodo di integrazione per la funzione level set di valore zero
    M_integrationMethod->set_integration_method ( mediumMesh.convex_index(), mediumIntegrationVelocity );
    M_integrationMethod->set_simplex_im ( integrationType );

    // Descrizione del metodo di integrazione per un lato della funzione level set
    M_integrationMethodInside->set_integration_method ( mediumMesh.convex_index(), mediumIntegrationVelocity );
    M_integrationMethodInside->set_simplex_im ( integrationType );

    // Descrizione del metodo di integrazione per l'altro lato della funzione level set
    M_integrationMethodOutside->set_integration_method ( mediumMesh.convex_index(), mediumIntegrationVelocity );
    M_integrationMethodOutside->set_simplex_im ( integrationType );

    M_levelSet->reinit();
    M_levelSet->values(1) = M_levelSet->values(0);


    // Riempio il level set con i valori nel file data
    for ( size_type d = 0; d < M_levelSet->get_mesh_fem().nb_basic_dof(); ++d )
    {
        const base_node node = M_levelSet->get_mesh_fem().point_of_basic_dof(d);

        M_levelSet->values(0) [ d ] = M_data->ylevelSetFunction(node);
        M_levelSet->values(1) [ d ] = M_data->levelSetCutFunction(node);
    }

    M_levelSet->touch();

    // Aggiorno la mesh e i metodi di integrazione poichè il level set è stato modificato
    M_mesh->adapt();
    M_integrationMethod->adapt();
    M_integrationMethodInside->adapt();
    M_integrationMethodOutside->adapt();

    // Valuto la funzione level set nei gradi di libertà primari
    const size_type shiftPressure = mediumMeshFEMPressure.nb_dof();
    // gives the total number of different degrees of freedom. If the optional reduction is used, this will be the number
	// of columns of the reduction matrix. Otherwise it will return the number of basic degrees of freedom

    gmm::resize(M_baricenterValue, shiftPressure);
    for ( size_type i = 0; i < shiftPressure; ++i )
    {
    	const base_node node = mediumMeshFEMPressure.point_of_basic_dof ( i );
    	M_baricenterValue [ i ] = M_data->ylevelSetFunction ( node );
    }

    // Valuto la funzione level set nei gradi di libertà secondari
    const size_type shiftVelocity = mediumMeshFEMVelocity.nb_dof();
    gmm::resize ( M_DOFValue, shiftVelocity );



    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
    	const base_node node = mediumMeshFEMVelocity.point_of_basic_dof ( i );
    	M_DOFValue [ i ] = M_data->ylevelSetFunction ( node );
    }

    return;
}// init
