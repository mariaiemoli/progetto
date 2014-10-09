/** LevelSetData.h
 *
 */

#ifndef LEVELSETDATA_H_
#define LEVELSETDATA_H_ 1

#include "Core.h"
#include "Parser.h"

class LevelSetData
{
public:

    LevelSetData ( const GetPot& dataFile,
                   const std::string& section = "fractureData/",
                   const std::string& sectionLevelSet = "levelSet/" );

    /** scalar_type ylevelSetFunction ( const base_node& x, int num = 0 )
     *
     * prima funzione level set che definisce la frattura, level set valutato in (x,y)
     *
     */
    scalar_type ylevelSetFunction ( const base_node& x,
                                   int num = 0 );


    /** scalar_type levelSetFunction ( const base_node& x, int num = 0 )
     *
     * seconda funzione level set che definisce la frattura, level set valutato in (t,y)
     *
     */
    scalar_type levelSetFunction ( const base_node& x,
                                   int num = 0 );


    /** scalar_type levelSetCutFunction ( const base_node& x, int num = 0 )
     *
     * terza funzione level set che definisce la frattura, se per caso voglio "tagliare la frattura"
     *
     */
    scalar_type levelSetCutFunction ( const base_node& x,
                                      int num = 0 );


    /** scalar_type y_map ( const base_node& t )
     *
     * mappa dalla frattura piatta a quella y(t)
     *
     */
    scalar_type y_map ( const base_node& t );


    /** scalar_type x_map ( const base_node& t )
     *
     * mappa dalla frattura piatta a quella x(t)
     *
     */
    scalar_type x_map ( const base_node& t );


    /** scalar_type t_map ( const base_node& t )
     *
     * mappa inversa per trovare t dato x
     *
     */
    scalar_type t_map ( const base_node& t );



    /** scalarVector_Type map_jac ( const base_node& x, const size_type& num );
     *
     * questa funzione serve per convertire la lunghezza/area degli elementi dalla frattura piatta a quella mappata
     *
     */
    scalarVector_Type map_jac ( const base_node& x,
                           const size_type& num );


    /**    scalarVector_Type normal_map ( const base_node& P, const size_type& num );
     *
     * questa funzione calcola la normale alla frattura - può anche dipendere da x
     *
     */
    scalarVector_Type normal_map ( const base_node& P,
                              const size_type& num );


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
    std::string M_map_inv;
    std::string M_map_jac;
    std::string M_normal_map;

    //SIMPLEX_INTEGRATION
    std::string M_integrationTypeSimplex;

    mutable LifeV::Parser M_parser;

};

typedef LevelSetData LevelSetData_Type;
typedef boost::shared_ptr<LevelSetData_Type> LevelSetDataPtr_Type;

#endif /* LEVELSETDATA_H_ */
