#include "../include/UsefulFunctions.h"

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
}

void exportSolution ( const std::string& fileName,
                      const std::string& solutionName,
                      const getfem::mesh_fem& meshFEM,
                      const scalarVector_Type& solution )
{
    getfem::vtk_export exporter(fileName.data());

    exporter.exporting(meshFEM);

    exporter.write_mesh();

    exporter.write_point_data(meshFEM, solution, solutionName.data());

}

void exportSolutionInCell ( const std::string& fileName,
                            const std::string& solutionName,
                            const getfem::mesh_fem& meshFEM,
                            const scalarVector_Type& solution )
{
    getfem::vtk_export exporter(fileName.data());

    exporter.exporting(meshFEM);

    exporter.write_mesh();

    exporter.write_cell_data(solution, solutionName.data());

}

void exportMesh ( const std::string& fileName, const getfem::mesh& mesh )
{
    getfem::vtk_export exporter(fileName.data());

    exporter.exporting(mesh);

    exporter.write_mesh();
}

scalar_type pointDistance ( const scalar_type& x0,
                            const scalar_type& x1,
                            const scalar_type& y0,
                            const scalar_type& y1 )
{
    return std::sqrt(std::pow(x0 - x1, 2) + std::pow(y0 - y1, 2));
}

void fromBitVectorToStdVector ( dal::bit_vector& bitVector, std::vector < size_type >& stdVector )
{
        size_type i_cv = 0;

        for ( i_cv << bitVector; i_cv != size_type(-1); i_cv << bitVector )
        {
                stdVector.push_back ( i_cv );
        }

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

bool isInTriangle ( const getfem::mesh& mesh, const size_type& elementID, const base_node& node, const scalar_type& toll )
{
    const base_node x0 = node - mesh.points_of_convex(elementID) [ 0 ];
    const base_node x1 = node - mesh.points_of_convex(elementID) [ 1 ];
    const base_node x2 = node - mesh.points_of_convex(elementID) [ 2 ];

    const scalar_type A = 0.5 * ( gmm::abs(x0 [ 0 ] * x1 [ 1 ] - x0 [ 1 ] * x1 [ 0 ]) +
                                  gmm::abs(x0 [ 0 ] * x2 [ 1 ] - x0 [ 1 ] * x2 [ 0 ]) +
                                  gmm::abs(x2 [ 0 ] * x1 [ 1 ] - x2 [ 1 ] * x1 [ 0 ]) );

    if ( gmm::abs( A - mesh.convex_area_estimate(elementID) ) <= toll )
    {
        return true;
    }
    else
    {
        return false;
    }

} // isInTriangle

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


std::pair < std::string, size_type >
comparaSegni ( const std::string& region, const scalarVector_Type& signs)
{
    std::pair < std::string, size_type > compara;
    compara.first = "Base";

    scalar_type flag = 0;

    for ( size_type i = 0; i < signs.size(); ++i )
    {
        if ( !( ( region[i] == '+' && signs[i] > 0 ) || ( region[i] == '-' && signs[i] < 0 ) ) )
        {
            ++flag;
            compara.first = "Extended";
            compara.second = i;
        }
    }

    if ( flag == signs.size() )
    {
        compara.first = "Extra";
    }

    return compara;

}

