#ifndef _FRACTUREINTERSECT_
#define _FRACTUREINTERSECT_ 1

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
                Cross = 500000
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

        size_type getNumberIntersectionOfType ( IntersectionType type ) const;

        size_type getNumberIntersections () const;

        size_type getNumberType () const;

        size_type getBasisFunctionOfType ( IntersectionType type ) const;

private:

        IntersectionType intersectionType ( getfem::mesh_level_set& meshLevelSet,
                                            stringContainer_Type& regionActive,
                                            const size_type& elementID,
                                            const sizeVector_Type& levelSets );

        scalar_type integrateWithBooleanOperation ( getfem::mesh_level_set& meshLevelSet,
                                                    const size_type& elementID,
                                                    const std::string& operation ) const;

        void findActiveRegion ( const scalarVector_Type& integrationValue,
                                stringContainer_Type& regionActive) const;

        mapIntersection_Type M_intersections;
        std::map < regionLevelSetPair_Type, IntersectionType > M_subRegionIntersection;
        std::map < IntersectionType, size_type> M_basisFunctionOfType;
        stringContainer_Type M_subRegion;

}; // class FractureIntersect

typedef FractureIntersect FractureIntersect_Type;
typedef boost::shared_ptr < FractureIntersect > FractureIntersectPtr_Type;

#endif // _FRACTUREINTERSECT_
