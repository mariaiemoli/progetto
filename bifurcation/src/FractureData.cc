

#include "../include/FractureData.h"


/**************************************************************************/
/*  FractureData.cc														  */
/*  Classe che contiene tutte le informazioni circa la natura geometrica  */
/*  e fisica della singola frattura										  */
/**************************************************************************/

FractureData::FractureData ( const GetPot& dataFile,
                             const std::string& section,
                             const std::string& sectionDomain,
                             const std::string& sectionDarcy ) :
            M_section ( section ),
            M_sectionDomain ( M_section + sectionDomain ),
            M_sectionDarcy ( M_section + sectionDarcy ),

            // domain
            M_thickness ( dataFile ( ( M_sectionDomain + "thickness" ).data (), 0.01 ) ),
            M_lengthAbscissa ( dataFile ( ( M_sectionDomain + "lengthAbscissa" ).data (), 1. ) ),
            M_lengthOrdinate ( dataFile ( ( M_sectionDomain + "lengthOrdinate" ).data (), 0. ) ),
            M_lengthQuota ( dataFile ( ( M_sectionDomain + "lengthQuota" ).data (), 0. ) ),
            M_spatialDiscretization ( dataFile ( ( M_sectionDomain + "spatialDiscretization" ).data (), 200 ) ),
            M_translateAbscissa ( dataFile ( ( M_sectionDomain + "translateAbscissa" ).data (), 0. ) ),
            M_spaceDimension ( dataFile ( ( M_section + "spaceDimension" ).data (), 1. ) ),
            M_integrationTypeVelocity ( dataFile ( ( M_sectionDomain + "integrationTypeVelocity" ).data (),
                                                   "IM_GAUSS1D(3)" ) ),
            M_integrationTypePressure ( dataFile ( ( M_sectionDomain + "integrationTypePressure" ).data (),
                                                   "IM_GAUSS1D(2)" ) ),
            M_meshType ( dataFile ( ( M_sectionDomain + "meshType" ).data (), "GT_PK(1,1)" ) ),
            M_fEMTypePressure ( dataFile ( ( M_sectionDomain + "FEMTypePressure" ).data (), "FEM_PK(1,0)" ) ),
            M_fEMTypeVelocity ( dataFile ( ( M_sectionDomain + "FEMTypeVelocity" ).data (), "FEM_PK(1,1)" ) ),
            M_fEMTypeLinear ( dataFile ( ( M_sectionDomain + "FEMTypeLinear" ).data (), "FEM_PK(1,1)" ) ),
            
            // darcy
            M_etaNormal ( dataFile ( ( M_sectionDarcy + "etaNormal" ).data (), M_thickness ) ),
            M_etaTangential ( dataFile ( ( M_sectionDarcy + "etaTangential" ).data (), 1. / M_thickness ) ),
            M_etaNormalDistribution ( dataFile ( ( M_sectionDarcy + "etaNormalDistribution" ).data (), "1." ) ),
            M_etaTangentialDistribution ( dataFile ( ( M_sectionDarcy + "etaTangentialDistribution" ).data (), "1." ) ),
            M_darcySource ( dataFile ( ( M_sectionDarcy + "source" ).data (), "1." ) ),
            M_pressureExact ( dataFile ( ( M_sectionDarcy + "solution" ).data (), "1." ) ),
            M_velocityExact ( dataFile ( ( M_sectionDarcy + "velocity" ).data (), "1." ) ),

            M_meshSpacing( dataFile ( ( M_sectionDomain + "spacing" ).data (), "x" ) )
{
}// costruttore


//questa bella funzione restituisce un'eventuale modulazione del coefficiente di permeabilità normale
scalar_type FractureData::etaNormalDistribution ( const base_node& x )
{
    M_parser.setString ( M_etaNormalDistribution );
    M_parser.setVariable ( "x", x [ 0 ] );
    
    return M_parser.evaluate ();
}// etaNormalDistribution

//questa bella funzione restituisce un'eventuale modulazione del coefficiente di permeabilità tangenziale
scalar_type FractureData::etaTangentialDistribution ( const base_node& x )
{
    M_parser.setString ( M_etaTangentialDistribution );
    M_parser.setVariable ( "x", x [ 0 ] );
    
    return M_parser.evaluate ();
}// etaTangentialDistribution


// Exact solution, velocity (non ho ancora impostato quella corretta)
scalar_type FractureData::velocityExact ( const base_node& x, const base_node& n )
{
    M_parser.setString ( M_velocityExact );
    M_parser.setVariable ( "x", x [ 0 ] );
    M_parser.setVariable ( "y", x [ 1 ] );
    M_parser.setVariable ( "n1", n [ 0 ] );
    M_parser.setVariable ( "n2", n [ 1 ] );
    
    return M_parser.evaluate ();
}// velocityExact

// Exact solution, div(Velocity) -- SET = 0 WITH NO MASS SOURCES/SINKS !
scalar_type FractureData::darcySource ( const base_node& x )
{
    M_parser.setString ( M_darcySource );
    M_parser.setVariable ( "x", x [ 0 ] );
    
    return M_parser.evaluate ();
}// darcySource

//pressione esatta nella frattura, se volessi fare il caso con p imposta - qui invece è predisposto per risolvere il problema accoppiato
scalar_type FractureData::pressureExact ( const base_node& x )
{
    M_parser.setString ( M_pressureExact );
    M_parser.setVariable ( "x", x [ 0 ] );
    
    return M_parser.evaluate ();
}// pressureExact


scalar_type FractureData::meshSpacing( const scalar_type& x )
{
    M_parser.setString ( M_meshSpacing );
    M_parser.setVariable ( "x", x );
    
    return M_parser.evaluate ();
}// meshSpacing
