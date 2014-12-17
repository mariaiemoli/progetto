
#include "../include/LevelSetData.h"

/**************************************************************************/
/*  LevelSetData.cc														  */
/*  Classe che contiene tutte le informazioni legate al level set e le 	  */
/*  funzioni per valutarne il valore nei punti                 			  */
/**************************************************************************/

LevelSetData::LevelSetData ( const GetPot& dataFile,
                             const std::string& section,
                             const std::string& sectionLevelSet ) :
            M_section ( section ),
            M_sectionLevelSet ( M_section + sectionLevelSet ),
            M_spaceDimension ( dataFile ( ( M_section + "spaceDimension" ).data (), 1. ) ),
            // level set
            M_function ( dataFile ( ( M_sectionLevelSet + "levelSet" ).data (), "x" ) ),
            M_xfunction ( dataFile ( ( M_sectionLevelSet + "xlevelSet" ).data (), "t" ) ),
            M_yfunction ( dataFile ( ( M_sectionLevelSet + "ylevelSet" ).data (), "t" ) ),
            M_cutFunction ( dataFile ( ( M_sectionLevelSet + "levelSetCut" ).data (), "-1" ) ),
            M_x_map ( dataFile ( ( M_sectionLevelSet + "xMap" ).data (), "1" ) ),
            M_y_map ( dataFile ( ( M_sectionLevelSet + "yMap" ).data (), "1" ) ),
            M_map_jac ( dataFile ( ( M_sectionLevelSet + "jacMap" ).data (), "1" ) ),
            M_normal_map ( dataFile ( ( M_sectionLevelSet + "normalMap" ).data (), "1" ) ),
            M_integrationTypeSimplex ( dataFile ( ( M_sectionLevelSet + "integrationTypeSimplex" ).data (),
                                                  "IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(3),1)" ) )
{}


scalar_type LevelSetData::ylevelSetFunction ( const base_node& t )
{
    M_parser.setString ( M_function );
    M_parser.setVariable ( "x", t [ 0 ] );
    M_parser.setVariable ( "y", t [ 1 ] );
    
    return M_parser.evaluate ();
}// ylevelSetFunction


scalar_type LevelSetData::levelSetFunction ( const base_node& t )
{
    M_parser.setString ( M_yfunction );
    M_parser.setVariable ( "t", t [ 0 ] );
    M_parser.setVariable ( "y", t [ 1 ] );
    
    return M_parser.evaluate ();
}// levelSetFunction


scalar_type LevelSetData::levelSetCutFunction ( const base_node& t )
{
    M_parser.setString ( M_cutFunction );
    M_parser.setVariable ( "t", t [ 0 ] );
    M_parser.setVariable ( "y", t [ 1 ] );
    
    return M_parser.evaluate ();
}// levelSetCutFunction


scalar_type LevelSetData::y_map ( const base_node& t )
{
    M_parser.setString ( M_y_map );
    M_parser.setVariable ( "t", t [ 0 ] );
    
    return M_parser.evaluate ();
}//  y_map


scalar_type LevelSetData::x_map ( const base_node& t )
{
    M_parser.setString ( M_x_map );
    M_parser.setVariable ( "t", t [ 0 ] );
    
    return M_parser.evaluate ();
}// x_map


scalarVector_Type LevelSetData::map_jac ( const base_node& x,
                                     const size_type& num )
{
    scalarVector_Type map_jac ( num, 0. );

    M_parser.setString ( M_map_jac );
    M_parser.setVariable ( "x", x [ 0 ] );
    M_parser.setVariable ( "y", x [ 1 ] );

    for ( size_type i = 0; i < num; ++i )
    {
        map_jac [ i ] = M_parser.evaluate ( i );
    }

    return map_jac;
}// map_jac


scalarVector_Type LevelSetData::normal_map ( const base_node& x,
                                        const size_type& num )
{
    scalarVector_Type normal_map ( num + 1, 0. );

    M_parser.setString ( M_map_jac );
    M_parser.setVariable ( "x", x [ 0 ] );
    M_parser.setVariable ( "y", x [ 1 ] );

    for ( size_type i = 0; i < num + 1; ++i )
    {
        normal_map [ i ] = M_parser.evaluate ( i );
    }

    return normal_map;
}// normal_map
