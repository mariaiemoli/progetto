/** FractureIntersect.h
 *
 */

#ifndef _FRACTUREINTERSECT_
#define _FRACTUREINTERSECT_ 1

#include <iostream>
#include "Core.h"
#include "IntersectData.h"
#include "FractureHandler.h"
#include <map>


class FractureIntersect
{
public:
        enum IntersectionType
        {
                Parallel = 400000,
                Cross = 500000,
                Bifurcation = 600000
        };


        typedef std::map < IntersectionType, IntersectDataContainer_Type > mapIntersection_Type;
        typedef std::pair < size_type, size_type > regionLevelSetPair_Type;


        FractureIntersect ();


        void constructIntesection ( getfem::mesh_level_set& meshLevelSet,
                                    const FracturePtrContainer_Type& fractures );

        IntersectDataContainer_Type& getIntersectionsOfType ( IntersectionType type )
        {
                return M_intersections [ type ];
        }

        mapIntersection_Type& getIntersections ()
        {
                return M_intersections;
        }


        /** size_type getNumberIntersectionOfType ( IntersectionType type ) const
         * restituisce il numero di intersezioni trovate del tipo type ( il numero di elementi tagliati)
         */
        size_type getNumberIntersectionOfType ( IntersectionType type ) const;
        
        size_type getNumberCross () const;
        
        size_type getNumberBifurcation () const;


        /** size_type getNumberIntersections () const
         * restituisce il numero totale di intersezioni
         */
        size_type getNumberIntersections () const;


        /**size_type getNumberType () const
         * restituisce il numero di intersezioni di tipo diverso trovate
         */
        size_type getNumberType () const;


        size_type getBasisFunctionOfType ( IntersectionType type ) const;



private:


        IntersectionType intersectionType ( getfem::mesh_level_set& meshLevelSet,
                                            const size_type& elementID,
                                            const sizeVector_Type& levelSets );


        scalar_type integrateWithBooleanOperation ( getfem::mesh_level_set& meshLevelSet,
                                                    const size_type& elementID,
                                                    const std::string& operation ) const;


        mapIntersection_Type M_intersections;
        std::map < regionLevelSetPair_Type, IntersectionType > M_subRegionIntersection;
		std::map < IntersectionType, size_type> M_basisFunctionOfType;

};

typedef FractureIntersect FractureIntersect_Type;
typedef boost::shared_ptr < FractureIntersect > FractureIntersectPtr_Type;

#endif /* _FRACTUREINTERSECT_H_ */
