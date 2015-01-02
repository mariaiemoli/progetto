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

#include "include/MeshHandler.h"
#include "include/FracturesSet.h"
#include "include/MediumData.h"
#include "include/DarcyFractured.h"


/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/


int main ( int argc, char* argv [ ] )
{
	std::string fileName("Bifurcation");
	
	if ( argc == 2 )
	{
		fileName = argv[1];
	}
	
	GetPot dataFile( fileName.c_str() );
	
	const std::string section = "";

	const std::string vtkFolder = "vtk/";

	// Creo la cartella dove salvare i risultati se già non esiste
	std::string s = "mkdir "+vtkFolder;
	system (s.c_str());

	//Data exporter
	std::cout << "Create the data exporter..." << std::flush;
	ExporterPtr_Type exporter( new Exporter_Type(dataFile));
	std::cout << " completed!" <<std::endl;


	//Mesh Handler
	std::cout << "Create the meshHandler..." << std::flush;
	MeshHandlerPtr_Type mesh(new MeshHandler_Type(dataFile, "mediumData/domain/"));
	mesh->setUpMesh();
	mesh->setUpFEM();
	std::cout<< " completed!" << std::endl;

	
	// Medium data for the Darcy problem
	std::cout << "Create the mediumData for the Darcy problem.." << std::flush;
	const std::string sectionSolverDarcy = "darcy/";
	MediumDataPtr_Type mediumDataDarcy(new MediumData_Type(dataFile,
										sectionSolverDarcy));
	std::cout << " completed!" << std::endl;
	

	// Fracture Set
	std::cout << "Create the set of fractures for " << std::flush;
	const size_type numberFractures = dataFile(
			(section + "numberFractures").data(), 0);
	std::cout << numberFractures << " fracture(s)..." << std::endl << std::flush;

	// Fracture Set
	FracturesSetPtr_Type fractures ( new FracturesSet );

	fractures->init ( dataFile, section, numberFractures, mesh->getMesh(), mesh->getMeshLevelSet(),
					  mesh->getIntegrationTypeVelocity(),
					  mesh->getMeshFEMScalar(), mesh->getMeshFEMVector() );

	std::cout << " completed!" << std::endl;


	// Create the mesh regions
	std::cout << "Create mesh regions..." << std::flush;
	mesh->setUpRegions ( fractures );
	std::cout << " completed!" << std::endl;

	// Fracture boundary conditions
	std::cout << "Create fracture boundary conditions..." << std::flush;
	BCPtrContainer_Type bcFracture(numberFractures);
	for ( size_type f = 0; f < numberFractures; ++f )
	{
		bcFracture [ f ].reset(new BC_Type( fractures->getFracture( f )->getMeshFlat(),
										   	fractures->getFracture ( f )->getData().getMeshType(),
										   	fractures->getFracture ( f ) -> getDofIntersection() ));
	}
	std::cout << " completed!" << std::endl;

	// Boundary conditions handler
	std::cout << "Create boundary conditions handler..." << std::flush;
	BCHandlerPtr_Type bcHandler(new BCHandler_Type( bcFracture ));
	std::cout << " completed!" << std::endl;


	// Save the cutted elements
	mesh->printCuttedElements(exporter->getFolder(), "cuttedElements.vtk");

	// Save the regions of the mesh
	exporter->meshRegion ( mesh->getMesh(), "RegionMesh.vtk" );
	

	// Compute inverse of mesh size (h^(-1) dove h è il passo di griglia)
	std::cout << "Compute inverse of mesh size..." << std::flush;
	mesh->computeMeshMeasures();
	for ( size_type f = 0; f < numberFractures; ++f )
	{
		fractures->getFracture( f )->computeInvH(bcHandler);
	}
	std::cout << " completed!" << std::endl;

	// Darcy problem
	std::cout << "Create Darcy problem..." << std::flush;
	DarcyFracturedPtr_Type darcy(new DarcyFractured_Type(mediumDataDarcy, mesh,
			bcHandler, fractures, exporter));
	std::cout << " completed!" << std::endl;

	// Initialize the solver
	std::cout << "Initialize the Darcy problem..." << std::flush;
	darcy->init();
	std::cout << " completed!" << std::endl;

	// Assembly the matrices and vectors
	std::cout << "Assembly the Darcy problem..." << std::flush;
	darcy->assembly( dataFile);
	std::cout << " completed!" << std::endl;

	// Solve and save the solutions
	std::cout << "Solve the Darcy problem..." << std::flush;
	darcy->solve();
	std::cout << " completed!" << std::endl;

	return 0;

}
