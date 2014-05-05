#include "../include/LevelSetHandler.h"

LevelSetHandler::LevelSetHandler ( const GetPot& dataFile,
                                   const std::string& section,
                                   const std::string& sectionLevelSet ) :
    M_data(new LevelSetData_Type(dataFile, section, sectionLevelSet))//, M_mesh(
//mesh), M_levelSet(mesh), M_data->getSpaceDimension(), true)

{
}

void LevelSetHandler::init ( getfem::mesh& mediumMesh,
                             const std::string& mediumIntegrationTypeVelocity,
                             const getfem::mesh_fem& mediumMeshFEMPressure,
                             const getfem::mesh_fem& mediumMeshFEMVelocity )
{
    M_mesh.reset(new GFMeshLevelSet_Type(mediumMesh));

    M_levelSet.reset(new GFLevelSet_Type(mediumMesh));

    M_integrationMethod.reset(new GFIntegrationMethodLevelSet_Type(*M_mesh,
            getfem::mesh_im_level_set::INTEGRATE_BOUNDARY));

    M_integrationMethodInside.reset(new GFIntegrationMethodLevelSet_Type(
            *M_mesh, getfem::mesh_im_level_set::INTEGRATE_INSIDE));

    M_integrationMethodOutside.reset(new GFIntegrationMethodLevelSet_Type(
            *M_mesh, getfem::mesh_im_level_set::INTEGRATE_OUTSIDE));

    // Level sets

    // Integration type for the level set function
    getfem::pintegration_method integrationType =
            getfem::int_method_descriptor(M_data->getIntegrationTypeSimplex());

    // Set that the M_mediumMesh M_levelSetMesh is cutted by a level set
    M_mesh->add_level_set(*M_levelSet);

    // Integration type for the dual variable
    getfem::pintegration_method mediumIntegrationVelocity =
            getfem::int_method_descriptor(mediumIntegrationTypeVelocity);

    // Set the integration method for the zero of the level set function
    M_integrationMethod->set_integration_method(mediumMesh.convex_index(),
            mediumIntegrationVelocity);

    M_integrationMethod->set_simplex_im(integrationType);

    // Set the integration method for one side of the level set function
    M_integrationMethodInside->set_integration_method(
            mediumMesh.convex_index(), mediumIntegrationVelocity);

    M_integrationMethodInside->set_simplex_im(integrationType);

    // Set the integration method for the other sife of the level set function
    M_integrationMethodOutside->set_integration_method(
            mediumMesh.convex_index(), mediumIntegrationVelocity);

    M_integrationMethodOutside->set_simplex_im(integrationType);

    M_levelSet->reinit();

    M_levelSet->values(1) = M_levelSet->values(0);
    // Fill the level set with the value from the data
    for ( size_type d = 0; d < M_levelSet->get_mesh_fem().nb_basic_dof(); ++d )
    {
        const base_node node = M_levelSet->get_mesh_fem().point_of_basic_dof(d);
        M_levelSet->values(0) [ d ] = M_data->levelSetFunction(node);
        M_levelSet->values(1) [ d ] = M_data->levelSetCutFunction(node);
    }

    M_levelSet->touch();

    // Since the level set is modified update the mesh and integration methods
    M_mesh->adapt();
    M_integrationMethod->adapt();
    M_integrationMethodInside->adapt();
    M_integrationMethodOutside->adapt();

    // Evaluate the level set function in the primal degrees of freedom
    const size_type shiftPressure = mediumMeshFEMPressure.nb_dof();
    gmm::resize(M_baricenterValue, shiftPressure);
    for ( size_type i = 0; i < shiftPressure; ++i )
    {
        M_baricenterValue [ i ] = M_data->levelSetFunction(
                mediumMeshFEMPressure.point_of_basic_dof(i), 0.0);
    }

    // Evaluate the level set function in the dual degrees of freedom
    const size_type shiftVelocity = mediumMeshFEMVelocity.nb_dof();
    gmm::resize(M_DOFValue, shiftVelocity);
    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        M_DOFValue [ i ] = M_data->levelSetFunction(
                mediumMeshFEMVelocity.point_of_basic_dof(i), 0.0);
    }

}
