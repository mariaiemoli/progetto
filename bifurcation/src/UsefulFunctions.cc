
#include "../include/UsefulFunctions.h"

/**************************************************************************/
/*  UsefulFunctions.cc													  */
/*  													                  */
/**************************************************************************/


void massLumping ( sparseMatrix_Type& matrix )
{

    const scalar_type matrixSize = gmm::mat_nrows(matrix);

    for ( size_type i = 0; i < matrixSize; ++i )
    {
        scalar_type lumped = 0;

        for ( size_type j = 0; j < matrixSize; ++j )
        {
            lumped += matrix(i, j);
            matrix(i, j) = 0;
        }
        matrix(i, i) = lumped;
    }
    
    return;
}// massLumping


void exportSolution ( const std::string& fileName,
                      const std::string& solutionName,
                      const getfem::mesh_fem& meshFEM,
                      const scalarVector_Type& solution )
{
    getfem::vtk_export exporter(fileName.data());

    exporter.exporting(meshFEM);

    exporter.write_mesh();

    exporter.write_point_data(meshFEM, solution, solutionName.data());
    
    return;

}// exportSolution


void exportSolutionInCell ( const std::string& fileName,
                            const std::string& solutionName,
                            const getfem::mesh_fem& meshFEM,
                            const scalarVector_Type& solution )
{
    getfem::vtk_export exporter(fileName.data());
 
    exporter.exporting(meshFEM);
 
    exporter.write_mesh();
 
    exporter.write_cell_data(solution, solutionName.data());
    
    return;
 
}// exportSolutionInCell


void exportMesh ( const std::string& fileName, const getfem::mesh& mesh )
{
    getfem::vtk_export exporter(fileName.data(), true);
    
    // con true mi permette di esportare la mesh in formato txt
     
    exporter.exporting(mesh);

    exporter.write_mesh();
    
    return;
}// exportMesh


scalar_type pointDistance ( const scalar_type& x0,
                            const scalar_type& x1,
                            const scalar_type& y0,
                            const scalar_type& y1 )
{
    return std::sqrt(std::pow(x0 - x1, 2) + std::pow(y0 - y1, 2));
}// pointDistance


void fromBitVectorToStdVector ( dal::bit_vector& bitVector, std::vector < size_type >& stdVector )
{
        size_type i_cv = 0;

        for ( i_cv << bitVector; i_cv != size_type(-1); i_cv << bitVector )
        {
                stdVector.push_back ( i_cv );
        }
        
        return;

} // fromBitVectorToStdVector


char intToChar ( const size_type& integer )
{
        return static_cast<char>( integer + 97 );
} // intToChar


std::string regionSigns ( const scalarVector_Type& levelSetValue )
{
    std::ostringstream regionSigns;
    const size_type numLevelSet = levelSetValue.size();
    
    for ( size_type i = 0; i < numLevelSet; ++i )
    {
        if ( levelSetValue[i] > 0 )
        {
            regionSigns << "+";
        }
        else
        {
            regionSigns << "-";
        }
    }

    return regionSigns.str();
} // regionSigns


std::string getOperation ( const std::string& subRegion, const sizeVector_Type& levelSets )
{
    std::ostringstream operation;

    for ( size_type i = 0; i < levelSets.size(); ++i )
    {
        if ( subRegion [i] == '-' )
        {
            operation << "(!" << intToChar ( levelSets[i] ) << ")";
        }
        else
        {
            operation << intToChar ( levelSets[i] );
        }

        // Not the last one
        if ( i != levelSets.size() - 1 )
        {
            operation << "*";
        }
    }

    return operation.str().c_str();
} // getOperation

void orderId( size_type& id_i, size_type& id_j, size_type& id_k )
{
	size_type id0, id1, id2;
	
	id0 = fmin( id_i, fmin( id_j, id_k ) );
	id2 = fmax( id_i, fmax( id_j, id_k ) );
	
	if( id0 == id_i )
	{
		id1 = fmin( id_j, id_k );
	}
	else if( id0 == id_j )
	{
		id1= fmin( id_i, id_k );
	}
	else
	{
		id1 = fmin( id_i, id_j );
	}
	
	id_i = id0;
	id_j = id1;
	id_k = id2;
	
	return;
	
} // orderId

