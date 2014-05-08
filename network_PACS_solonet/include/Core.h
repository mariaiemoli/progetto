/*! Questo file:
	 - include tutte le librerie getfem che ci servono
	 - definisce tutti i typedef che poi useremo
	 - definisce le enumeraizone utilizzate più avanti nel codice
*/

#ifndef _DARCYCORE_
#define _DARCYCORE_ 1

// Level Set and Xfem stuff:
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_interpolation.h>
#include <getfem/getfem_config.h>
#include <getfem/getfem_assembling.h> // import assembly methods (and comp. of 
#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>     // export functions (save the solution in a file)
#include <gmm/gmm.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_MUMPS_interface.h>
#include <gmm/gmm_superlu_interface.h>

#include <getfem/bgeot_mesh.h>

#include <getfem/getfem_mesh_im_level_set.h>
#include <getfem/getfem_mesh_fem_level_set.h>

#include <string>

#include <getfem/getfem_mesh_fem_product.h>
#include <getfem/getfem_mesh_fem_global_function.h>

#include "GetPot"
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <iomanip>

/*! bgeot:: è un namespace di una libreria di  getfem
*/

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector;
/* special class for small (dim < 16) vectors */
using bgeot::base_node;
/* geometrical nodes (derived from base_small_vector)*/
using bgeot::scalar_type;
/* = double */
using bgeot::size_type;
/* = unsigned long */
//using bgeot::short_type;

/* definition of some matrix/vector types. These ones are built
 using the predefined types in Gmm++ */

typedef gmm::rsvector<scalar_type> sparseVector_Type;

typedef gmm::row_matrix<sparseVector_Type> sparseMatrix_Type;
typedef boost::shared_ptr<sparseMatrix_Type> sparseMatrixPtr_Type;
typedef std::vector<sparseMatrix_Type> sparseMatrixContainer_Type;
typedef std::vector<sparseMatrixPtr_Type> sparseMatrixPtrContainer_Type;

typedef std::vector<scalar_type> scalarVector_Type;
typedef boost::shared_ptr<scalarVector_Type> scalarVectorPtr_Type;
typedef std::vector<scalarVector_Type> scalarVectorContainer_Type;
typedef boost::shared_ptr<scalarVectorPtr_Type> scalarVectorContainerPtr_Type;
typedef std::vector<scalarVectorPtr_Type> scalarVectorPtrContainer_Type;

typedef std::vector<size_type> sizeVector_Type;
typedef boost::shared_ptr<sizeVector_Type> sizeVectorPtr_Type;
typedef std::vector<sizeVector_Type> sizeVectorContainer_Type;
typedef std::vector<sizeVectorPtr_Type> sizeVectorPtrContainer_Type;

typedef std::pair<size_type, size_type> pairSize_Type;
typedef std::vector<pairSize_Type > pairSizeVector_Type;
typedef std::vector<pairSizeVector_Type> pairSizeVectorContainer_Type;

typedef std::vector < std::string > stringContainer_Type;

typedef getfem::mesh_fem GFMeshFEM_Type;
typedef boost::shared_ptr<GFMeshFEM_Type> GFMeshFEMPtr_Type;
typedef std::vector<GFMeshFEMPtr_Type> GFMeshFEMPtrContainer_Type;

typedef getfem::mesh_level_set GFMeshLevelSet_Type;
typedef boost::shared_ptr<GFMeshLevelSet_Type> GFMeshLevelSetPtr_Type;
typedef std::vector < GFMeshLevelSetPtr_Type > GFMeshLevelSetPtrContainer_Type;

typedef getfem::level_set GFLevelSet_Type;
typedef std::vector < GFLevelSet_Type > GFLevelSetContainer_Type;
typedef boost::shared_ptr<GFLevelSet_Type> GFLevelSetPtr_Type;
typedef std::vector < GFLevelSetPtr_Type > GFLevelSetPtrContainer_Type;

typedef getfem::mesh_im_level_set GFIntegrationMethodLevelSet_Type;
typedef boost::shared_ptr<GFIntegrationMethodLevelSet_Type> GFIntegrationMethodLevelSetPtr_Type;

namespace LifeV
{
typedef scalar_type Real;
typedef size_type UInt;
typedef size_type ID;
typedef int Int;
}

enum ElementDimension
{
    MEDIUM = 2, FRACTURE = MEDIUM-1
};

/*! introduco anche i typedef contenuti in Parser messi in un commento in modo che possiamo trovarli già nel Core.h senza doverci ricordare le dipendenza che ci portano a Parser.h
 *
 * Parser:  
 *	typedef std::vector< std::string >                       stringsVector_Type;
 *  typedef std::string::const_iterator                      stringIterator_Type;
 *  typedef ParserSpiritGrammar< stringIterator_Type >       calculator_Type;
 *  typedef calculator_Type::results_Type                    results_Type;
 *  typedef IteratorType                                        iterator_Type;
 *  typedef boost::iterator_range< iterator_Type >              iteratorRange_Type;
 *  typedef ResultsType                                         results_Type;
 *
 * BC
 * typedef BC BC_Type;
 * typedef boost::shared_ptr<BC> BCPtr_Type;
 * typedef std::vector<BCPtr_Type> BCPtrContainer_Type;
 *
 * LevelSetHandler
 * typedef LevelSetHandler LevelSetHandler_Type;
 * typedef boost::shared_ptr<LevelSetHandler_Type> LevelSetHandlerPtr_Type;
 *
 * FractureHandler
 * typedef FractureHandler FractureHandler_Type;
 * typedef boost::shared_ptr<FractureHandler_Type> FractureHandlerPtr_Type;
 * typedef std::vector<FractureHandlerPtr_Type> FracturePtrContainer_Type;
 * typedef std::vector<FractureHandler_Type> FractureContainer_Type;
 * typedef boost::shared_ptr<FracturePtrContainer_Type> FracturePtrContainerPtr_Type;
 *
 *FractureIntersect
 * enum IntersectionType
 *      {
 *              Parallel = 400000,
 *              Cross = 500000
 *      };
 * typedef std::map < IntersectionType, IntersectDataContainer_Type > mapIntersection_Type;
 * typedef std::pair < size_type, size_type > regionLevelSetPair_Type;
 *
*/

#endif
