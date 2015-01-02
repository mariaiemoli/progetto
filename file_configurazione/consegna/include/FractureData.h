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


#ifndef _FRACTUREDATA_
#define _FRACTUREDATA_ 1

#include "Core.h"
#include "Parser.h"

/**************************************************************************/
/*  FractureData.h														  */
/*  Classe che contiene tutte le informazioni circa la natura geometrica  */
/*  e fisica della singola frattura										  */
/**************************************************************************/


class FractureData
{
public:
    enum
    {
        FRACTURE = 0
    };

    FractureData ( const GetPot& dataFile,
                   const std::string& section = "fractureData/",
                   const std::string& sectionDomain = "domain/",
                   const std::string& sectionDarcy = "darcy/" );


    /**
     * Funzione che restituisce un'eventuale modulazione del coefficiente di permeabilità in direzione normale.
     */
    scalar_type etaNormalDistribution ( const base_node& x );

    /**
     * Funzione che restituisce un'eventuale modulazione del coefficiente di permeabilità in direzione tangenziale.
     */
    scalar_type etaTangentialDistribution ( const base_node& x );

    // Soluzione esatta, velocità 
    scalar_type velocityExact ( const base_node& x,
                                const base_node& n );

   
    // Soluzione esatta, div(Velocity) -- SET = 0 WITH NO MASS SOURCES/SINKS !
    scalar_type darcySource ( const base_node& x );

    scalar_type pressureExact ( const base_node& x );

    scalar_type meshSpacing ( const scalar_type& x );

    inline scalar_type getThickness () const
    {
        return M_thickness;
    }


    inline scalar_type getEtaNormal () const
    {
        return M_etaNormal;
    }


    inline scalar_type getEtaTangential () const
    {
        return M_etaTangential;
    }


    inline scalar_type getLengthAbscissa () const
	{
		return M_lengthAbscissa;
	}


	inline scalar_type getLengthOrdinate () const
	{
		return M_lengthOrdinate;
	}


	inline scalar_type getLengthQuota () const
	{
		return M_lengthQuota;
	}


    inline size_type getSpatialDiscretization () const
    {
        return M_spatialDiscretization;
    }


    inline scalar_type getTranslateAbscissa () const
    {
        return M_translateAbscissa;
    }


    inline bgeot::dim_type getSpaceDimension () const
    {
        return M_spaceDimension;
    }


    inline std::string getIntegrationTypeVelocity () const
    {
        return M_integrationTypeVelocity;
    }


    inline std::string getIntegrationTypePressure () const
    {
        return M_integrationTypePressure;
    }


    inline std::string getMeshType () const
    {
        return M_meshType;
    }


    inline std::string getFEMTypePressure () const
    {
        return M_fEMTypePressure;
    }


    inline std::string getFEMTypeVelocity () const
    {
        return M_fEMTypeVelocity;
    }


    inline std::string getFEMTypeLinear () const
    {
        return M_fEMTypeLinear;
    }


private:

    // Attributes
    std::string M_section;
    std::string M_sectionDomain;
    std::string M_sectionDarcy;

    // Proprietà del mezzo

    scalar_type M_thickness; // thickness of the fracture

    scalar_type M_lengthAbscissa;
	scalar_type M_lengthOrdinate;
	scalar_type M_lengthQuota;

    size_type M_spatialDiscretization;

    scalar_type M_translateAbscissa;

    bgeot::dim_type M_spaceDimension;

    std::string M_integrationTypeVelocity;
    std::string M_integrationTypePressure;

    std::string M_meshType;

    std::string M_fEMTypePressure;
    std::string M_fEMTypeVelocity;
    std::string M_fEMTypeLinear;

    scalar_type M_etaNormal; // eta_gamma = d/K_normale
    scalar_type M_etaTangential; // eta_t=1/(K_t*d)
    std::string M_etaNormalDistribution;
    std::string M_etaTangentialDistribution;

    std::string M_darcySource;
    std::string M_pressureExact;
    std::string M_velocityExact;

    std::string M_meshSpacing;

    mutable LifeV::Parser M_parser;

};

#endif /* FRACTUREDATA_H_ */
