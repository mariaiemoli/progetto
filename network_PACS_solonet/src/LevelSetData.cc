#include "../include/LevelSetData.h"

LevelSetData::LevelSetData ( const GetPot& dataFile,
                             const std::string& section,
                             const std::string& sectionLevelSet ) :
            M_section ( section ),												// numero frattura
            M_sectionLevelSet ( M_section + sectionLevelSet ),
            M_spaceDimension ( dataFile ( ( M_section + "spaceDimension" ).data (), 1. ) ),	// se non è definito assegno valore 1. di default
            // level set
            M_function ( dataFile ( ( M_sectionLevelSet + "levelSet" ).data (), "x" ) ),
            M_cutFunction ( dataFile ( ( M_sectionLevelSet + "levelSetCut" ).data (), "-1" ) ),
            M_z_map ( dataFile ( ( M_sectionLevelSet + "zMap" ).data (), "1" ) ),
            M_map_jac ( dataFile ( ( M_sectionLevelSet + "jacMap" ).data (), "1" ) ),
            M_normal_map ( dataFile ( ( M_sectionLevelSet + "normalMap" ).data (), "1" ) ),
            M_integrationTypeSimplex ( dataFile ( ( M_sectionLevelSet + "integrationTypeSimplex" ).data (),
                                                  "IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(3),1)" ) )
{
}
/*
 * M_section 			-->  è una stringa, rappresenta semplicemente il nome della frattura nel file data
 * M_sectionLevelSet	-->  è una stringa, rappresenta il nome della sezione relativa alla descrizione della funzione level set
 * M_spaceDimention		-->  intero di dimensione fissa
 * M_function			-->  è una stringa, gli passo il campo levelSet del file data, se non è definita assegna di default valore x
 * M_cutFunction		-->  stringa
 * M_z_map				-->  stringa, rappresenta la mappa per la variabile z
 * M_map_jac			-->  stringa
 * M_normal_map			-->  stringa
 * M_integrationTypeSimplex	-->  stringa
 *
 * l'unico campo che non viene riempito dal costruttore è il campo M_parser, sono le funzioni che definiscono la frattura
 * che si occupano di riempirlo
 */

// prima funzione level set che definisce la frattura . decommentare quello voluto
scalar_type LevelSetData::levelSetFunction ( const base_node& x,
                                             int num )
{
    M_parser.setString ( M_function );		// definisco la funzione che rappresenta la frattura
    M_parser.setVariable ( "x", x [ 0 ] );
    M_parser.setVariable ( "y", x [ 1 ] );
    return M_parser.evaluate ();
}

/* LifeV::Parser M_parser
 * oggetto della classe Parser
 */

//seconda funzione level set (di finestratura come dico io :D) se per caso voglio "tagliare la frattura"
scalar_type LevelSetData::levelSetCutFunction ( const base_node& x,
                                                int num )
{
    M_parser.setString ( M_cutFunction );
    M_parser.setVariable ( "x", x [ 0 ] );
    M_parser.setVariable ( "y", x [ 1 ] );
    return M_parser.evaluate ();
}

//mappa dalla frattura piatta a quella y(x) - decommentare il caso voluto
scalar_type LevelSetData::z_map ( const base_node& x )
{
    M_parser.setString ( M_z_map );
    M_parser.setVariable ( "x", x [ 0 ] );
    return M_parser.evaluate ();
}

// questo serve per convertire la lunghezza/area degli elementi dalla frattura piatta a quella mappata
// decommentare il caso voluto. può anche dipendere da x volendo
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
}

//normale alla frattura - può anche dipendere da x
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
}

