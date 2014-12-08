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


#ifndef USEFULFUNCTIONS_H_
#define USEFULFUNCTIONS_H_ 1

#include "Core.h"

/**************************************************************************/
/*  UsefulFunctions.cc													  */
/*  													                  */
/**************************************************************************/


/** 
 * Funzione che approssima una matrice sparsa in una matrice diagonale non singolare con la tecnica mass lumping che permette di 
 * disaccoppiare tra loro le equazioni del sistema
 */
void massLumping ( sparseMatrix_Type& matrix );


/**
 * Funzioni per esportare soluzioni e mesh
 */
void exportSolution ( const std::string& fileName,
                      const std::string& solutionName,
                      const getfem::mesh_fem& meshFEM,
                      const scalarVector_Type& solution );

void exportSolutionInCell ( const std::string& fileName,
                            const std::string& solutionName,
                            const getfem::mesh_fem& meshFEM,
                            const scalarVector_Type& solution );

void exportMesh ( const std::string& fileName, const getfem::mesh& mesh );


/**
 * Funzione che calcola la distanza tra due punti
 */
scalar_type pointDistance ( const scalar_type& x0,
                            const scalar_type& x1,
                            const scalar_type& y0,
                            const scalar_type& y1 );


void fromBitVectorToStdVector ( dal::bit_vector& bitVector, std::vector < size_type >& stdVector );


char intToChar ( const size_type& integer );


/**
 * Funzione che restituisce ' + ' se il levelset ha valore positivo, ' - ' in caso contrario
 */
std::string regionSigns ( const scalarVector_Type& levelSetValue );


/**
 * Funzione che restituisce un'operazione tra due levelset
 */
std::string getOperation ( const std::string& subRegion, const sizeVector_Type& levelSets );

void orderId( size_type& id_i, size_type& id_j, size_type& id_k );

					

#endif /* USEFULFUNCTIONS_H_ */
