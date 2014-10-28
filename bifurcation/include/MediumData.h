/**
 *
 * MediumData.h
 *
 * classe che contiene le informazioni sulle proprietà del mezzo e sui metodi di integrazione
 * 
 */

#ifndef _MEDIUMDATA_
#define _MEDIUMDATA_ 1

#include "Core.h"
#include "Parser.h"

class MediumData
{
public:

	/**
	 * \param GetPot& dataFile: variabile di tipo GetPot che contiene il nome del file data
	 * \param std::string& section: nome della sezione corrispondente al mezzo
	 * \param std::string& section = "mediumData/": nome della sezione del mezzo corrispondente alle caratteristiche del dominio,
	 * 												se non fornito è posto di default pari a mediumData/"
	 */
	MediumData ( const GetPot& dataFile,
				 const std::string& sectionSolver,
				 const std::string& section = "mediumData/",
				 const std::string& sectionDomain = "domain/" );


    /** 
     * Funzione che valuta la soluzione esatta in un dato nodo.
     * \param base_node& x: bgeot::base_node, nodo geometrico
     *
     */
    scalar_type exact ( const base_node& x, const scalar_type& t = 0 ) const;


    /** 
     * Funzione che valuta la soluzione esatta in un dato nodo nella regione in cui level set > 0. 
     * \param base_node& x: bgeot::base_node, nodo geometrico
     */
    scalar_type exactOutlet ( const base_node& x, const scalar_type& t = 0 ) const;



    /** 
     * Funzione che valuta la soluzione esatta in un dato nodo nella regione in cui level set < 0. 
     * \param base_node& x: bgeot::base_node, nodo geometrico
     */
    scalar_type exactInlet ( const base_node& x, const scalar_type& t = 0 ) const;


    /** 
     * Funzione che valuta il flusso esatto in un dato nodo.  
     * \param base_node& x: bgeot::base_node, nodo geometrico
     * \param base_node& n: bgeot::base_node, vettore normale
     */
    scalar_type exactFlux ( const base_node& x,
                            const base_node& n,
                            const scalar_type& t = 0 ) const;


    /** 
     * Funzione che valuta il termine sorgente in un dato nodo. 
     * \param base_node& x: bgeot::base_node, nodo geometrico
     */
    scalar_type source ( const base_node& x, const scalar_type& t = 0 ) const;


    /**
     * Funzione che restituisce un'eventuale modulazione del coefficiente di permeabilità.
     * \param base_node& x: bgeot::base_node, nodo geometrico
     */
    scalar_type invKDistribution11 ( const base_node& x ) const;


    /**
     * Funzione che restituisce un'eventuale modulazione del coefficiente di permeabilità.
     * \param base_node& x: bgeot::base_node, nodo geometrico
     */
    scalar_type invKDistribution12 ( const base_node& x ) const;


    /**
     * Funzione che restituisce un'eventuale modulazione del coefficiente di permeabilità.
     * \param base_node& x: bgeot::base_node, nodo geometrico
     */
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

typedef MediumData MediumData_Type;									/*!< classe MediumData */
typedef boost::shared_ptr<MediumData_Type> MediumDataPtr_Type;		/*!< puntatore alla classe MediumData */ 


#endif /* MEDIUMDATA_H_ */
