#ifndef LEVELSETHANDLER_H_
#define LEVELSETHANDLER_H_ 1

#include "LevelSetData.h"

class LevelSetHandler
{
public:

    LevelSetHandler ( const GetPot& dataFile, const std::string& section =
            "fractureData/", const std::string& sectionLevelSet = "levelSet/" );

    void
    init ( getfem::mesh& mediumMesh,
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
    // The level sets defining the fracture
    GFLevelSetPtr_Type M_levelSet;

    //valore del level set nei baricentri
    scalarVector_Type M_baricenterValue;
    //valore del level set nei gradi di libertà di v
    scalarVector_Type M_DOFValue;

    // integration method, level set (on the boundary)
    GFIntegrationMethodLevelSetPtr_Type M_integrationMethod;
    // integration method, level set (inside)
    GFIntegrationMethodLevelSetPtr_Type M_integrationMethodInside;
    // integration method, level set (outside)
    GFIntegrationMethodLevelSetPtr_Type M_integrationMethodOutside;

};

typedef LevelSetHandler LevelSetHandler_Type;
typedef boost::shared_ptr<LevelSetHandler_Type> LevelSetHandlerPtr_Type;

#endif /* LEVELSETHANDLER_H_ */
