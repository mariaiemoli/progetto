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


#ifndef _FRACTURESSET_
#define _FRACTURESSET_ 1

#include "Core.h"
#include "FractureHandler.h"
#include "FractureIntersect.h"

/**************************************************************************/
/*  FractureSet.h														  */
/*  Classe che contiene tutte le fratture				                  */
/**************************************************************************/


class FracturesSet
{
public:
	
	// Costruttore nullo
	FracturesSet ();


	/**
	 * Funzione che inizializza l'insieme delle fratture definendo un vettore con le fratture e una classe delle intersezioni.
	 */
    void init ( const GetPot& dataFile,
                const std::string& section,
                const size_type& numFractures,
                getfem::mesh& mesh,
                getfem::mesh_level_set& meshLevelSet,
                const std::string& integrationTypeVelocity,
                const getfem::mesh_fem& meshFEMScalar,
                const getfem::mesh_fem& meshFEMVector );


    /**
     * Funzione che restituisce il numero totale di fratture.
     */
    size_type getNumberFractures () const
    {
            return M_fractures.size();
    }


    /**
     * Funzione che restituisce il puntatore ad una frattura.
     * \param size_type& f: id della frattura che mi interessa
     * \return FractureHandlerPtr_Type& M_fractures[f]: frattura di indice " f " nell'insieme delle fratture
     */
    const FractureHandlerPtr_Type& getFracture ( const size_type& f ) const
    {
            return M_fractures[f];
    }


    /**
     * Funzione che restituisce il puntatore ad una frattura.
     * \param size_type& f: id della frattura che mi interessa
     * \return FractureHandlerPtr_Type& M_fractures[f]: frattura di indice " f " nell'insieme delle fratture
     */
    FractureHandlerPtr_Type& getFracture ( const size_type& f )
    {
            return M_fractures[f];
    }


    /**
     * Funzione che retituisce l'insieme di tutte le intersezioni.
     * \return FractureIntersectPtr_Type& M_intersections: puntatore alla classe delle intersezioni
     */
    const FractureIntersectPtr_Type& getIntersections () const
    {
            return M_intersections;
    }


    /**
     * Funzione che retituisce l'insieme di tutte le intersezioni.
     * \return FractureIntersectPtr_Type& M_intersections: puntatore alla classe delle intersezioni
     */

    FractureIntersectPtr_Type& getIntersections ()
    {
            return M_intersections;
    }


private:
	FracturePtrContainer_Type M_fractures;

	FractureIntersectPtr_Type M_intersections;

};


typedef FracturesSet FracturesSet_Type;										/*!< Classe FracturesSet */
typedef boost::shared_ptr < FracturesSet_Type > FracturesSetPtr_Type;		/*!< Puntatore alla classe FracturesSet */


#endif /* _FRACTURESSET_H_ */
