/**
 * Interior Penalty for Darcy's equation
 *
 */
#include "include/MeshHandler.h"
#include "include/FracturesSet.h"
#include "include/MediumData.h"
#include "include/DarcyFractured.h"
#include "include/TransportFractured.h"

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main ( int argc, char* argv [ ] )
{

    GetPot command_line(argc, argv);

    const std::string data_file_name = command_line.follow("data", 2, "-f",
            "--file");

    GetPot dataFile(data_file_name.data());

    const std::string section = "";

    const std::string vtkFolder = "vtk/";

    // Data exporter
    std::cout << "Create the data exporter..." << std::flush;
    ExporterPtr_Type exporter(new Exporter_Type(dataFile));
    std::cout << " completed!" << std::endl;

    // Mesh handler
    std::cout << "Create the meshHandler..." << std::flush;
    MeshHandlerPtr_Type mesh(new MeshHandler_Type(dataFile,
            "mediumData/domain/"));
    mesh->setUpMesh();
    mesh->setUpFEM();
    std::cout << " completed!" << std::endl;

    // Medium data for the Darcy problem
    std::cout << "Create the mediumData for the Darcy problem..." << std::flush;
    const std::string sectionSolverDarcy = "darcy/";
    MediumDataPtr_Type mediumDataDarcy(new MediumData_Type(dataFile,
            sectionSolverDarcy));
    std::cout << " completed!" << std::endl;

    // Fracture Set
    std::cout << "Create the set of fractures for " << std::flush;
    const size_type numberFractures = dataFile(
            (section + "numberFractures").data(), 0);
    std::cout << numberFractures << " fracture(s)..." << std::flush;

    // Create the set of the fractures
    FracturesSetPtr_Type fractures ( new FracturesSet );

    fractures->init ( dataFile, section, numberFractures, mesh->getMesh(), mesh->getMeshLevelSet(),
                      mesh->getIntegrationTypeVelocity(),
                      mesh->getMeshFEMScalar(), mesh->getMeshFEMVector() );

    std::cout << " completed!" << std::endl;

    // Create the mesh regions
    std::cout << "Create mesh regions..." << std::flush;
    mesh->setUpRegions ( fractures );
    std::cout << " completed!" << std::endl;

    // Medium boundary conditions
    std::cout << "Create medium boundary conditions..." << std::flush;
    BCPtr_Type bcMedium(new BC_Type(mesh->getMesh(), mesh->getMeshType(),
            MEDIUM));
    std::cout << " completed!" << std::endl;

    // Fracture boundary conditions
    std::cout << "Create fracture boundary conditions..." << std::flush;
    BCPtrContainer_Type bcFracture(numberFractures);
    for ( size_type f = 0; f < numberFractures; ++f )
    {
        bcFracture [ f ].reset(new BC_Type(fractures->getFracture( f )->getMeshFlat(),
                                           fractures->getFracture ( f )->getData().getMeshType(),
                                           FRACTURE));
    }
    std::cout << " completed!" << std::endl;

    // Boundary conditions handler
    std::cout << "Create boundary conditions handler..." << std::flush;
    BCHandlerPtr_Type bcHandler(new BCHandler_Type(bcMedium, bcFracture));
    std::cout << " completed!" << std::endl;

    std::cout << "Setup boundary conditions handler..." << std::flush;
    bcHandler->createBDRegions(mesh->getMesh());
    bcHandler->createBDRegionsFractures(mesh->getMesh());
    std::cout << " completed!" << std::endl;

    // Save the cutted elements
    mesh->printCuttedElements(exporter->getFolder(), "cuttedElements.vtk");

    // Save the regions of the mesh
    exporter->meshRegion ( mesh->getMesh(), "RegionMesh.vtk" );

    // Compute inverse of mesh size
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
    darcy->assembly();
    std::cout << " completed!" << std::endl;

    // Solve and save the solutions
    std::cout << "Solve the Darcy problem..." << std::flush;
    darcy->solve();
    std::cout << " completed!" << std::endl;

  
    return 0;
}
