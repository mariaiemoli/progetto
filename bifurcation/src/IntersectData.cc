/** IntersectData.cc
 *
 */

#include "../include/IntersectData.h"

void IntersectData::copy ( const IntersectData& in )
{
    if ( this != &in )
    {
        M_fractures.resize ( in.M_fractures.size() );
        for ( size_type i = 0; i < M_fractures.size(); ++i )
        {
            M_fractures[i] = in.M_fractures[i];
        }

        M_elementID = in.M_elementID;

        M_dofPressure.resize ( in.M_dofPressure.size() );
        for ( size_type i = 0; i < M_dofPressure.size(); ++i )
        {
            M_dofPressure[i] = in.M_dofPressure[i];
        }

        M_dofVelocity.resize ( in.M_dofVelocity.size() );
        for ( size_type i = 0; i < M_dofVelocity.size(); ++i )
        {
            M_dofVelocity[i] = in.M_dofVelocity[i];
        }

        M_matrices = in.M_matrices;
    }
} // copy


void IntersectData::setIntersection ( const size_type& elementID,
									  const FracturePtrContainer_Type& fractures, const size_type k)
{
    M_elementID = elementID;
    M_fractures = fractures;
    
    if ( k == 1) 	// k==1 -> Bifurcation
    {
    	M_matrices.setMatrices ( M_fractures );
    }
    
} // setIntersection

