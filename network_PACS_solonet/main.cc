/**
 * Interior Penalty for Darcy's equation
 *
 */
#include "include/MeshHandler.h"			// cartella include
#include "include/FracturesSet.h"			// cartella include
#include "include/MediumData.h"				// cartella include
#include "include/DarcyFractured.h"			// cartella include
#include "include/TransportFractured.h"		// cartella include

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main ( int argc, char* argv [ ] )
{
	/*
	 * leggo i dati in input usando getpot
	 *I parametri vengono letti mediante l’operatore (), associato ad un oggetto
	 *di tipo GetPot, in cui speciﬁchiamo il nome dato al parametro e
	 *l’eventuale valore di default.
	 *./main a=10. b=70. nint=100  <- esempio di chiamata da terminale
	 *argc -> numero parametri in ingresso
	 *argv -> lista parametri in ingresso
	 */

	// uso getpot sia con argc e argv che con la lettura da file

    GetPot command_line(argc, argv);

    // nel file getpot-doc nella cartella dropbox è spiegato il significato
    // se trova -f o --file sostituisce quello che segue con data
    const std::string data_file_name = command_line.follow("data", 2, "-f",
            "--file");

    GetPot dataFile(data_file_name.data());	// oggetto di tipo getpot letto da file

    const std::string section = "";

    const std::string vtkFolder = "vtk/";

    //////////////////////////////////////////////////////////
    // Data exporter		--------> Exporter.h
    std::cout << "Create the data exporter..." << std::flush;

    /*
     * flush() makes sure the file stream is updated with the data.
	 * endl() puts a newline and then uses flush().
	 * flush è uno stream di output
     */

    // definito in Exporter.h, puntatore alla classe Exporter
    ExporterPtr_Type exporter(new Exporter_Type(dataFile));
    std::cout << " completed!" << std::endl;


    //////////////////////////////////////////////////////////
    // Mesh handler			--------> MeshHander.h
    std::cout << "Create the meshHandler..." << std::flush;
    // definito in MeshHandler.h, puntatore alla classe MeshHandler
    MeshHandlerPtr_Type mesh(new MeshHandler_Type(dataFile,
            "mediumData/domain/"));
    mesh->setUpMesh();
    mesh->setUpFEM();
    std::cout << " completed!" << std::endl;

    //////////////////////////////////////////////////////////
    // Medium data for the Darcy problem		---------> MediumData.h
    std::cout << "Create the mediumData for the Darcy problem..." << std::flush;
    const std::string sectionSolverDarcy = "darcy/";
    // definito in MediumData.h, puntatore alla classe MediumData
    MediumDataPtr_Type mediumDataDarcy(new MediumData_Type(dataFile,
            sectionSolverDarcy));
    std::cout << " completed!" << std::endl;

    //////////////////////////////////////////////////////////
    // Fracture Set		--------> leggo da file con un oggetto di tipo getpot
    std::cout << "Create the set of fractures for " << std::flush;
    const size_type numberFractures = dataFile(
            (section + "numberFractures").data(), 0);
    std::cout << numberFractures << " fracture(s)..." << std::flush;

    //////////////////////////////////////////////////////////
    // Create the set of the fractures		 -------> FracturesSet.h
    FracturesSetPtr_Type fractures ( new FracturesSet );

    fractures->init ( dataFile, section, numberFractures, mesh->getMesh(), mesh->getMeshLevelSet(),
                      mesh->getIntegrationTypeVelocity(),
                      mesh->getMeshFEMScalar(), mesh->getMeshFEMVector() );

    std::cout << " completed!" << std::endl;


    //////////////////////////////////////////////////////////
    // Create the mesh regions		------> MeshHandler.h
    std::cout << "Create mesh regions..." << std::flush;
    mesh->setUpRegions ( fractures );
    std::cout << " completed!" << std::endl;


    //////////////////////////////////////////////////////////
    // Medium boundary conditions		--------> BC.h
    std::cout << "Create medium boundary conditions..." << std::flush;
    BCPtr_Type bcMedium(new BC_Type(mesh->getMesh(), mesh->getMeshType(),
            MEDIUM));
    std::cout << " completed!" << std::endl;


    //////////////////////////////////////////////////////////
    // Fracture boundary conditions		-------> BC.h
    std::cout << "Create fracture boundary conditions..." << std::flush;
    BCPtrContainer_Type bcFracture(numberFractures);
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        bcFracture [ f ].reset(new BC_Type(fractures->getFracture( f )->getMeshFlat(),
                                           fractures->getFracture ( f )->getData().getMeshType(),
                                           FRACTURE));
    }
    std::cout << " completed!" << std::endl;


    //////////////////////////////////////////////////////////
    // Boundary conditions handler			-------> BCHandler.h
    std::cout << "Create boundary conditions handler..." << std::flush;
    BCHandlerPtr_Type bcHandler(new BCHandler_Type(bcMedium, bcFracture));
    std::cout << " completed!" << std::endl;

    std::cout << "Setup boundary conditions handler..." << std::flush;
    bcHandler->createBDRegions(mesh->getMesh());
    bcHandler->createBDRegionsFractures(mesh->getMesh());
    std::cout << " completed!" << std::endl;


    //////////////////////////////////////////////////////////
    // Save the cutted elements			--------> MeshHandler.h
    mesh->printCuttedElements(exporter->getFolder(), "cuttedElements.vtk");


    //////////////////////////////////////////////////////////
    // Save the regions of the mesh		---------> Exporter.h
    exporter->meshRegion ( mesh->getMesh(), "RegionMesh.vtk" );


    //////////////////////////////////////////////////////////
    // Compute inverse of mesh size			--------> MeshHandler.h, FractureHandler.h
    std::cout << "Compute inverse of mesh size..." << std::flush;
    mesh->computeMeshMeasures();
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        fractures->getFracture( f )->computeInvH(bcHandler);
    }
    std::cout << " completed!" << std::endl;


    //////////////////////////////////////////////////////////
    // Darcy problem			-------> DarcyFractured.h
    std::cout << "Create Darcy problem..." << std::flush;
    DarcyFracturedPtr_Type darcy(new DarcyFractured_Type(mediumDataDarcy, mesh,
            bcHandler, fractures, exporter));
    std::cout << " completed!" << std::endl;


    //////////////////////////////////////////////////////////
    // Initialize the solver		-------> DarcyFractured.h
    std::cout << "Initialize the Darcy problem..." << std::flush;
    darcy->init();
    std::cout << " completed!" << std::endl;


    //////////////////////////////////////////////////////////
    // Assembly the matrices and vectors	-------> DarcyFractured.h
    std::cout << "Assembly the Darcy problem..." << std::flush;
    darcy->assembly();
    std::cout << " completed!" << std::endl;


    //////////////////////////////////////////////////////////
    // Solve and save the solutions		-------> DarcyFractured.h
    std::cout << "Solve the Darcy problem..." << std::flush;
    darcy->solve();
    std::cout << " completed!" << std::endl;

  
    return 0;
}
