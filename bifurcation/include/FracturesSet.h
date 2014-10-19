/**
 *
 * FractureSet.h
 * 
 * classe che contiene tutte le fratture
 *
 */

#ifndef _FRACTURESSET_
#define _FRACTURESSET_ 1

#include "Core.h"
#include "FractureHandler.h"
#include "FractureIntersect.h"


class FracturesSet
{
public:
	
	FracturesSet ();


	/**
	 * funzione che inizializza l'insieme delle fratture definendo un vettore con le fratture e una classe delle intersezioni
	 */
    void init ( const GetPot& dataFile,
                const std::string& section,
                const size_type& numFractures,
                getfem::mesh& mesh,
                getfem::mesh_level_set& meshLevelSet,
                const std::string& integrationTypeVelocity,
                const getfem::mesh_fem& meshFEMScalar,
                const getfem::mesh_fem& meshFEMVector );


    size_type getNumberFractures () const
    {
            return M_fractures.size();
    }


    const FractureHandlerPtr_Type& getFracture ( const size_type& f ) const
    {
            return M_fractures[f];
    }


    FractureHandlerPtr_Type& getFracture ( const size_type& f )
    {
            return M_fractures[f];
    }


    const FractureIntersectPtr_Type& getIntersections () const
    {
            return M_intersections;
    }


    FractureIntersectPtr_Type& getIntersections ()
    {
            return M_intersections;
    }


private:
	FracturePtrContainer_Type M_fractures;

	FractureIntersectPtr_Type M_intersections;

};


typedef FracturesSet FracturesSet_Type;
typedef boost::shared_ptr < FracturesSet_Type > FracturesSetPtr_Type;


#endif /* _FRACTURESSET_H_ */
