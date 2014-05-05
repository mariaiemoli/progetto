#ifndef _MEDIUMDATA_
#define _MEDIUMDATA_ 1

#include "Core.h"
#include "Parser.h"

class MediumData
{
public:

    MediumData ( const GetPot& dataFile,
                 const std::string& sectionSolver,
                 const std::string& section = "mediumData/",
                 const std::string& sectionDomain = "domain/" );

    // Exact solution, pressure
    scalar_type exact ( const base_node& x, const scalar_type& t = 0 ) const;

    // Exact solution, pressure OUT, i.e. level set >0
    scalar_type
    exactOutlet ( const base_node& x, const scalar_type& t = 0 ) const;

    // Exact solution, pressure IN (i.e. level set <0)
    scalar_type
    exactInlet ( const base_node& x, const scalar_type& t = 0 ) const;

    // Exact solution, velocity (non ho ancora impostato quella corretta)
    scalar_type exactFlux ( const base_node& x,
                            const base_node& n,
                            const scalar_type& t = 0 ) const;

    // Exact solution, div(Velocity) -- SET = 0 WITH NO MASS SOURCES/SINKS !
    scalar_type source ( const base_node& x, const scalar_type& t = 0 ) const;

    //questa bella funzione restituisce un'eventuale modulazione del coefficiente di permeabilità
    scalar_type invKDistribution11 ( const base_node& x ) const;

    scalar_type invKDistribution12 ( const base_node& x ) const;

    scalar_type invKDistribution22 ( const base_node& x ) const;

    // Initial solution for time dependent problem
    scalar_type exactInitial ( const base_node& x ) const;

    inline scalar_type getInvK ( ) const
    {
        return M_invK;
    }

    inline scalar_type getPenaltyVector ( ) const
    {
        return M_penaltyVector;
    }

    inline scalar_type getPenaltyScalar ( ) const
    {
        return M_penaltyScalar;
    }

private:

    std::string M_section;
    std::string M_sectionDomain;
    std::string M_sectionSolver;

    // Penalty parameters
    scalar_type M_penaltyVector;
    scalar_type M_penaltyScalar;

    // inverse permeability of the medium
    scalar_type M_invK;
    std::string M_invKDistribution11;
    std::string M_invKDistribution12;
    std::string M_invKDistribution22;
    std::string M_exact;
    std::string M_exactInlet;
    std::string M_exactOutlet;
    std::string M_exactFlux;
    std::string M_source;
    std::string M_exactInitial;

    mutable LifeV::Parser M_parser;

};

typedef MediumData MediumData_Type;
typedef boost::shared_ptr<MediumData_Type> MediumDataPtr_Type;

#endif
