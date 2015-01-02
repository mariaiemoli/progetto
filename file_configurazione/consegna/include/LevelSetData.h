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


#ifndef LEVELSETDATA_H_
#define LEVELSETDATA_H_ 1

#include "Core.h"
#include "Parser.h"

/**************************************************************************/
/*  LevelSetData.h														  */
/*  Classe che contiene tutte le informazioni legate al level set e le 	  */
/*  funzioni per valutarne il valore nei punti                 			  */
/**************************************************************************/


class LevelSetData
{
public:

	// Costruttore
    LevelSetData ( const GetPot& dataFile,
                   const std::string& section = "fractureData/",
                   const std::string& sectionLevelSet = "levelSet/" );

    /** 
     *
     * Questa funzione definisce il level set che rappresenta la frattura, level set valutato in (x,y).
     * \param base_node& x: nodo in coordinate ( x, y ) in cui valutare il levelset
     * \return scalar_type: valore del levelset nel nodo x
     *
     */
    scalar_type ylevelSetFunction ( const base_node& x );


    /** 
     *
     * Questa funzione definisce il level set che rappresenta la frattura, level set valutato in (t,y).
     * \param base_node& x: nodo in coordinate ( t, y ) in cui valutare il levelset
     * \return scalar_type: valore del levelset nel nodo x
     */
    scalar_type levelSetFunction ( const base_node& x );


    /** scalar_type levelSetCutFunction ( const base_node& x, int num = 0 )
     *
     * Questa funzione definisce il  level set che rappresenta la frattura, se per caso voglio "tagliare la frattura".
     * \param base_node& x: nodo in coordinate ( t, y ) in cui valutare il levelset
     * \return scalar_type: valore del levelset nel nodo x
     */
    scalar_type levelSetCutFunction ( const base_node& x );


    /** 
     * Questa funzione rappresenta la mappa dalla frattura piatta, y(t).
     */
    scalar_type y_map ( const base_node& t );


    /** 
	 * Questa funzione rappresenta la mappa dalla frattura piatta, x(t).
	 */
    scalar_type x_map ( const base_node& t );


    /** 
     * Questa funzione serve per convertire la lunghezza/area degli elementi dalla frattura piatta a quella mappata.
     */
    scalarVector_Type map_jac ( const base_node& x, const size_type& num );


    /** 
     * Funzione che calcola la normale alla frattura - può anche dipendere da x
     */
    scalarVector_Type normal_map ( const base_node& P, const size_type& num );


    inline bgeot::dim_type getSpaceDimension () const
    {
        return M_spaceDimension;
    }


    inline std::string getIntegrationTypeSimplex () const
    {
        return M_integrationTypeSimplex;
    }


    std::string getLevelSetFunctionString () const
    {
        return M_yfunction;
    }

private:

    // Attributes
    std::string M_section;
    std::string M_sectionLevelSet;

    bgeot::dim_type M_spaceDimension;

    std::string M_function;
    std::string M_xfunction;
    std::string M_yfunction;
    std::string M_cutFunction;

    std::string M_x_map;
    std::string M_y_map;
    std::string M_map_jac;
    std::string M_normal_map;

    //SIMPLEX_INTEGRATION
    std::string M_integrationTypeSimplex;

    mutable LifeV::Parser M_parser;

};

typedef LevelSetData LevelSetData_Type;								/*!< Classe LevelSetData */
typedef boost::shared_ptr<LevelSetData_Type> LevelSetDataPtr_Type;	/*!< Puntatore alla classe LevelSetData */

#endif /* LEVELSETDATA_H_ */
