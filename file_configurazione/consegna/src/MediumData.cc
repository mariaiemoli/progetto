
#include "../include/MediumData.h"

/**************************************************************************/
/*  MediumData.cc											              */
/*  Classe che contiene le informazioni sulle proprietà del mezzo e       */
/*  sui metodi di integrazione                 							  */
/**************************************************************************/

MediumData::MediumData ( const GetPot& dataFile,
                         const std::string& sectionSolver,
                         const std::string& section,
                         const std::string& sectionDomain ) :
            M_section( section ),
            M_sectionDomain( M_section + sectionDomain ),
            M_sectionSolver( M_section + sectionSolver ),
            M_penaltyVector( dataFile ( ( M_sectionDomain + "penaltyVelocity" ).data(), 5. ) ),
            M_penaltyScalar( dataFile ( ( M_sectionDomain + "penaltyPressure" ).data(), 1. ) ),
            // darcy 
            M_invK( dataFile ( ( M_sectionSolver + "invK" ).data(), 1. ) ),
            M_invKDistribution11( dataFile ( ( M_sectionSolver + "invKDist11" ).data(), "1." ) ),
            M_invKDistribution12( dataFile ( ( M_sectionSolver + "invKDist12" ).data(), "1." ) ),
            M_invKDistribution22( dataFile ( ( M_sectionSolver + "invKDist22" ).data(), "1." ) )
{
}// costruttore


scalar_type MediumData::invKDistribution11 ( const base_node& x ) const
{
    M_parser.setString(M_invKDistribution11);

    M_parser.setVariable("x", x [ 0 ]);
    M_parser.setVariable("y", x [ 1 ]);

    return M_parser.evaluate();
}// invKDistribution11


//questa bella funzione restituisce un'eventuale modulazione del coefficiente di permeabilità
scalar_type MediumData::invKDistribution12 ( const base_node& x ) const
{
    M_parser.setString(M_invKDistribution12);

    M_parser.setVariable("x", x [ 0 ]);
    M_parser.setVariable("y", x [ 1 ]);

    return M_parser.evaluate();
}// invKDistribution12


//questa bella funzione restituisce un'eventuale modulazione del coefficiente di permeabilità
scalar_type MediumData::invKDistribution22 ( const base_node& x ) const
{
    M_parser.setString(M_invKDistribution22);

    M_parser.setVariable("x", x [ 0 ]);
    M_parser.setVariable("y", x [ 1 ]);

    return M_parser.evaluate();
}// invKDistribution22

