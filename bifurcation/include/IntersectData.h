/** IntersectData.h
 *
 */

#ifndef _INTERSECTDATA_
#define _INTERSECTDATA_ 1

#include "FractureHandler.h"

class IntersectData
{

public:

    IntersectData ()
    {
        ;
    } // costruttore nullo


    IntersectData ( const IntersectData& in )
    {
        copy ( in );
    } // costruttore di copia


    void setIntersection ( const size_type& elementID,
                           const FracturePtrContainer_Type& fractures)
    {
        M_elementID = elementID;
        M_fractures = fractures;
    } // setIntersection


    void setDOFPosition ( const sizeVector_Type& dofPressure,
                          const sizeVectorContainer_Type& dofVelocity )
    {
        M_dofPressure = dofPressure;
        M_dofVelocity = dofVelocity;
    } // setDofPosition


    size_type getDOFPressure ( const size_type& component ) const
    {
        return M_dofPressure [ component ];
    } // getDOFPressure


    size_type getDOFVelocity ( const size_type& component, const size_type& region ) const
    {
        return M_dofVelocity [ component ] [ region ];
    } //  getDOFVelocity


    const sizeVector_Type& getDOFVelocity ( const size_type& region ) const
    {
        return M_dofVelocity [ region ];
    } // getDOFVelocity


    const FracturePtrContainer_Type& getFractures () const
    {
        return M_fractures;
    } // getFractures


    const FractureHandlerPtr_Type& getFracture ( const size_type& f ) const
    {
        return M_fractures [ f ];
    } // getFracture


    size_type getNumFractures () const
    {
        return M_fractures.size();
    } // getNumFractures


    const size_type& getElementID () const
    {
        return M_elementID;
    } // getElementID


    const IntersectData& operator = ( const IntersectData & in )
    {
        copy ( in );
    } // operator =


private:

    void copy ( const IntersectData& in );

    FracturePtrContainer_Type M_fractures;
    size_type M_elementID;
    sizeVector_Type M_dofPressure;
    sizeVectorContainer_Type M_dofVelocity;

};

typedef IntersectData IntersectData_Type;
typedef std::vector < IntersectData_Type > IntersectDataContainer_Type;
typedef boost::shared_ptr < IntersectData_Type > IntersectDataPtr_Type;
typedef std::vector < IntersectDataPtr_Type > IntersectDataPtrContainer_Type;


#endif /* _INTERSECTDATA_H_ */
