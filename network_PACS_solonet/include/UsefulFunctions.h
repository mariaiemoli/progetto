#ifndef USEFULFUNCTIONS_H_
#define USEFULFUNCTIONS_H_ 1

#include "Core.h"

void exportSolution ( const std::string& fileName,
                      const std::string& solutionName,
                      const getfem::mesh_fem& meshFEM,
                      const scalarVector_Type& solution );

void exportSolutionInCell ( const std::string& fileName,
                            const std::string& solutionName,
                            const getfem::mesh_fem& meshFEM,
                            const scalarVector_Type& solution );

void exportMesh ( const std::string& fileName, const getfem::mesh& mesh );

void massLumping ( sparseMatrix_Type& matrix );

void fromBitVectorToStdVector ( dal::bit_vector& bitVector,
                                std::vector < size_type >& stdVector );

char intToChar ( const size_type& integer );

std::string regionSigns ( const scalarVector_Type& levelSetValue );

scalar_type pointDistance ( const scalar_type& x0,
                            const scalar_type& x1,
                            const scalar_type& y0,
                            const scalar_type& y1 );

std::string getOperation ( const std::string& subRegion,
                           const sizeVector_Type& levelSets );

std::pair < std::string, size_type > comparaSegni ( const std::string& region,
                                                    const scalarVector_Type& signs);

bool isInTriangle ( const getfem::mesh& mesh,
                    const size_type& elementID,
                    const base_node& node,
                    const scalar_type& toll = 1e-7 );

#endif /* USEFULFUNCTIONS_H_ */
