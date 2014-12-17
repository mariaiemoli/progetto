
#include "../include/Exporter.h"

/**************************************************************************/
/*  Exporter.cc															  */
/*  Classe in cui definisco le funzioni per esportare i dati              */
/**************************************************************************/


Exporter::Exporter ( const GetPot& dataFile, const std::string& section ) :
					M_vtkFolder(dataFile((section + "folderVTK").data(), "./vtk/"))
{
}// costruttore


void Exporter::spy ( const sparseMatrixPtr_Type& matrix, const std::string& nameFile ) const
{
    gmm::MatrixMarket_IO::write(nameFile.c_str(), *matrix);

    return;
    
} // spy


void Exporter::spy ( const scalarVectorPtr_Type& vector, const std::string& nameFile ) const
{
    std::ofstream file;
    file.open(nameFile.c_str());
    for ( size_type i = 0; i < vector->size(); ++i )
    {
        file << std::setprecision(15) << vector->at(i) << std::endl;
    }
    file.close();
    
    return;
    
} // spy


void Exporter::meshRegion ( const getfem::mesh& mesh, const std::string& nameFile ) const
{
	// imposto lo spazio degli elementi finiti P0, ogni elemento contiene il valore della region mesh 
	GFMeshFEM_Type meshFEM ( mesh );
	getfem::pfem FEType = getfem::fem_descriptor ( "FEM_PK(2,0)" );
	meshFEM.set_finite_element ( FEType );

	// vettore che contiene i valori della region mesh, un valore per ogni grado di libertà
	const size_type nbElements = meshFEM.nb_basic_dof();
	scalarVector_Type regionMesh ( nbElements, 0 );

	// estraggo i flags delle rispettive regioni dalla mesh
	dal::bit_vector regionsIndexBitVector ( mesh.regions_index() );
	sizeVector_Type regionsIndex;
	fromBitVectorToStdVector ( regionsIndexBitVector, regionsIndex );

	for ( size_type i = 0; i < regionsIndex.size(); ++i )
	{
		getfem::mesh_region meshCurrentRegion ( mesh.region ( regionsIndex [i] ) );
		
		// per ogni regione della mesh controllo se è costituita solo da convessi, cioè se non contiene bordi di nessun convesso
		if ( meshCurrentRegion.is_only_convexes() == true )
		{
			dal::bit_vector currentRegionBitVector ( meshCurrentRegion.index() );
			size_type i_cv = 0;
			for ( i_cv << currentRegionBitVector; i_cv != bgeot::size_type(-1); i_cv << currentRegionBitVector )
			{
				const size_type position = (meshFEM.ind_basic_dof_of_element( i_cv ))[0];
				regionMesh [ position ] += regionsIndex [i];
			}
		}
	}

	// esporto la soluzione
	exportSolution ( M_vtkFolder + nameFile, "RegionMesh", meshFEM, regionMesh );
	
	return;

} // meshRegion
