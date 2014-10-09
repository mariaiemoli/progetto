#ifndef _INTERSECTDATA_
#define _INTERSECTDATA_ 1

#include "FractureHandler.h"

class IntersectData
{

public:

    IntersectData ()
    {
        ;
    }

    IntersectData ( const IntersectData& in )
    {
        copy ( in );
    }

    void setIntersection ( const size_type& elementID,
                           const FracturePtrContainer_Type& fractures,
                           const stringContainer_Type& regionActive,
                           const stringContainer_Type& subRegion )
    {
        M_elementID = elementID;
        M_fractures = fractures;
        M_regionActive = regionActive;
        findFacingRegion ( subRegion );
    }

    void setDOFPosition ( const sizeVector_Type& dofPressure,
                          const sizeVectorContainer_Type& dofVelocity )
    {
        M_dofPressure = dofPressure;
        M_dofVelocity = dofVelocity;
    }

    size_type getDOFPressure ( const size_type& component ) const
    {
        return M_dofPressure [ component ];
    }

    size_type getDOFVelocity ( const size_type& component, const size_type& region ) const
    {
        return M_dofVelocity [ component ] [ region ];
    }

    const sizeVector_Type& getDOFVelocity ( const size_type& region ) const
    {
        return M_dofVelocity [ region ];
    }

    const FracturePtrContainer_Type& getFractures () const
    {
        return M_fractures;
    }

    const FractureHandlerPtr_Type& getFracture ( const size_type& f ) const
    {
        return M_fractures [ f ];
    }

    size_type getNumFractures () const
    {
        return M_fractures.size();
    }

    const size_type& getElementID () const
    {
        return M_elementID;
    }

    const stringContainer_Type& getRegionActive () const
    {
        return M_regionActive;
    }

    size_type getIndexRegion ( const std::string& region ) const
    {
        return size_type (  std::find( M_regionActive.begin(), M_regionActive.end(), region ) - M_regionActive.begin() );
    }

    void findFacingRegion ( const stringContainer_Type& subRegion );

    const sizeVector_Type& getFacingRegion ( const size_type& i ) const
    {
        return M_facingRegion[i];
    }

    const IntersectData& operator = ( const IntersectData & in )
    {
        copy ( in );
    }

private:

    void copy ( const IntersectData& in );

    FracturePtrContainer_Type M_fractures;
    size_type M_elementID;
    stringContainer_Type M_regionActive;
    sizeVector_Type M_dofPressure;
    sizeVectorContainer_Type M_dofVelocity;
    sizeVectorContainer_Type M_facingRegion;


}; // class IntersectData

typedef IntersectData IntersectData_Type;
typedef std::vector < IntersectData_Type > IntersectDataContainer_Type;
typedef boost::shared_ptr < IntersectData_Type > IntersectDataPtr_Type;
typedef std::vector < IntersectDataPtr_Type > IntersectDataPtrContainer_Type;

#endif // _INTERSECTDATA_
