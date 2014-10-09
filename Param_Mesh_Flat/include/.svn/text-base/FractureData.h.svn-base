#ifndef _FRACTUREDATA_
#define _FRACTUREDATA_ 1

#include "Core.h"
#include "Parser.h"

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
                   const std::string& sectionDarcy = "darcy/",
                   const std::string& sectionTransport = "transport/" );

    // termine sorgente per il tracciante, problema di trasporto
    scalar_type transportSource ( const base_node& x,
                                  const scalar_type& time );

    // condizioni iniziali per il problema della concentrazione
    scalar_type concentrationInitialCondition ( const base_node& x );

    //questa bells funzione restituisce un'eventuale modulazione del coefficiente 
    //di permeabilità in direzione normale
    scalar_type etaNormalDistribution ( const base_node& x );

    //questa bells funzione restituisce un'eventuale modulazione del coefficiente 
    //di permeabilità in direzione tangenziale
    scalar_type etaTangentialDistribution ( const base_node& x );

    scalar_type muNormalDistribution ( const base_node& x );

    scalar_type muTangentialDistribution ( const base_node& x );

    // Exact solution, concentration - dipende anche dal tempo!
    scalar_type concentrationExact ( const base_node& x,
                                     const scalar_type& time );

    // Exact solution, concentration IN/OUT con una flag
    scalar_type Concentration_inout ( const base_node& x,
                                      const scalar_type& time,
                                      const size_type flag );

    // Exact solution, velocity (non ho ancora impostato quella corretta)
    scalar_type velocityExact ( const base_node& x,
                                const base_node& n );

    // Exact solution, flusso per la concentrazione - dipende anche dal tempo!
    scalar_type fluxExact ( const base_node& x,
                            const base_node& n,
                            const scalar_type& time );

    // Exact solution, div(Velocity) -- SET = 0 WITH NO MASS SOURCES/SINKS !
    scalar_type darcySource ( const base_node& x );

    scalar_type pressureExact ( const base_node& x );

    scalar_type meshSpacing ( const scalar_type& x );

    inline scalar_type getPosition () const
    {
        return M_position;
    }

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

    inline scalar_type getMuNormal () const
    {
        return M_muNormal;
    }

    inline scalar_type getMuTangential () const
    {
        return M_muTangential;
    }

    inline scalar_type getCsi0 () const
    {
        return M_csi0;
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
    std::string M_sectionTransport;

    scalar_type M_position; //posizione della frattura nel caso di frattura orizzontale a y=M_fracturePosition

    //properties of the media

    scalar_type M_thickness; // thickness of the fracture

    scalar_type M_csi0; //csi-0.5

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

    //inverse diffusivities
    scalar_type M_muNormal; // eta_gamma = d/K_normale
    scalar_type M_muTangential; // eta_t=1/(K_t*d)
    std::string M_muNormalDistribution;
    std::string M_muTangentialDistribution;

    std::string M_transportSource;
    std::string M_concentrationInitialCondition;

    std::string M_concentrationExact;
    std::string M_Concentration_inout;
    std::string M_fluxExact;

    std::string M_meshSpacing;

    mutable LifeV::Parser M_parser;

};

#endif
