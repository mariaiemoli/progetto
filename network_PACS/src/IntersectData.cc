#include "../include/IntersectData.h"

void IntersectData::findFacingRegion ( const stringContainer_Type& subRegion )
{
    std::vector < stringContainer_Type > facingRegion;

    M_facingRegion.resize(2);
    facingRegion.resize(2);

    M_facingRegion[0].resize ( 2 * ( M_regionActive.size() - 2 ) );
    M_facingRegion[1].resize ( 2 * ( M_regionActive.size() - 2 ) );

    facingRegion[0].resize ( 2 * ( M_regionActive.size() - 2 ) );
    facingRegion[1].resize ( 2 * ( M_regionActive.size() - 2 ) );

    sizeVector_Type pairFound (2, 0);

    for ( size_type i = 0; i < M_regionActive.size(); ++i )
    {
        const std::string first = M_regionActive[i];
        std::string::const_iterator it_first;

        for ( size_type j = i + 1; j < M_regionActive.size(); ++j )
        {
            const std::string second = M_regionActive[j];
            std::string::const_iterator it_second = second.begin();
            sizeVector_Type isDifferent(2,0);
            size_type count = 0;

            for ( it_first = first.begin(); it_first != first.end(); ++it_first, ++it_second, ++count )
            {
                if ( *it_first != *it_second )
                {
                    isDifferent [ count ] = 1;
                }
            }

            if ( isDifferent[0] + isDifferent[1] == 1 )
            {
                size_type indexFracture = 0;

                if ( isDifferent[1] == 1 )
                {
                    indexFracture = 1;
                }

                if ( first[indexFracture] == '-' )
                {
                    facingRegion [ indexFracture ] [ 2*pairFound[indexFracture] ] = first;
                    facingRegion [ indexFracture ] [ 2*pairFound[indexFracture] + 1 ] = second;
                }
                else
                {
                    facingRegion [ indexFracture ] [ 2*pairFound[indexFracture] ] = second;
                    facingRegion [ indexFracture ] [ 2*pairFound[indexFracture] + 1 ] = first;
                }

                (pairFound[indexFracture])++;

            }

        }
    }

    if ( facingRegion[0].size() > 2 )
    {
        for ( size_type indexFracture = 0; indexFracture < 2; ++indexFracture )
        {
            if ( facingRegion[indexFracture][0][(indexFracture+1)%2] != '-')
            {
                std::swap ( facingRegion[indexFracture][0], facingRegion[indexFracture][2]);
                std::swap ( facingRegion[indexFracture][1], facingRegion[indexFracture][3]);
            }
        }
    }

    for ( size_type i = 0; i < 2; ++i )
    {
        for ( size_type j = 0; j < facingRegion[i].size(); ++j )
        {
            //M_facingRegion[i][j] = size_type (std::find ( subRegion.begin(), subRegion.end(), facingRegion[i][j] ) - subRegion.begin());
            M_facingRegion[i][j] = size_type (std::find ( M_regionActive.begin(), M_regionActive.end(), facingRegion[i][j] ) - M_regionActive.begin());
        }
    }

} // findFacingRegion

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

        M_regionActive.resize( in.M_regionActive.size() );
        for ( size_type i = 0; i < M_regionActive.size(); ++i )
        {
            M_regionActive[i] = in.M_regionActive[i];
        }

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

        M_facingRegion.resize ( in.M_facingRegion.size() );
        for ( size_type i = 0; i < M_facingRegion.size(); ++i )
        {
            M_facingRegion [i].resize ( in.M_facingRegion[i].size() );
            for ( size_type j = 0; j < M_facingRegion[i].size(); ++j )
            {
                M_facingRegion[i][j] = in.M_facingRegion[i][j];
            }
        }
    }
} // copy

