/*
 * TransportFractured.h
 *
 *  Created on: Apr 1, 2011
 *      Author: fumagalli
 */

#ifndef TRANSPORTFRACTURED_H_
#define TRANSPORTFRACTURED_H_ 1

#include "Core.h"

#include "FractureHandler.h"
#include "XFEMOperators.h"
#include "BCHandler.h"
#include "MediumData.h"
#include "Exporter.h"
#include "TimeData.h"
#include "TransportStabilization.h"
#include "DarcyFractured.h"

class TransportFractured
{

public:
    TransportFractured ( const MediumDataPtr_Type& medium,
                         const MeshHandlerPtr_Type& mesh,
                         const BCHandlerPtr_Type& bcHandler,
                         FracturePtrContainer_Type& fracture,
                         const TimeDataPtr_Type& time,
                         const ExporterPtr_Type& exporter,
                         const DarcyFracturedPtr_Type& darcy );

    void init ( );

    void assembly ( );

    void addStabilization( const TransportStabilizationPtr_Type& stabilization );

    void solve (  );

private:

    // Data medium
    MediumDataPtr_Type M_mediumData;

    // Time data
    TimeDataPtr_Type M_timeData;

    // Mesh
    MeshHandlerPtr_Type M_mesh;

    // BC Handler
    BCHandlerPtr_Type M_bcHandler;

    // Fractures
    FracturePtrContainer_Type& M_fracture;

    // Exporter
    ExporterPtr_Type M_exporter;

    // Darcy solver
    DarcyFracturedPtr_Type M_darcy;

    //inverse diffusivities
    scalarVectorContainer_Type M_mediumMuInterpolated;
    // eta_gamma = d/K_normale - vector - sui punti della M_mediumMesh grande
    scalarVectorContainer_Type M_fractureMuNormalOnMedium;

    // Global matrix, advection-diffusion
    sparseMatrixPtr_Type M_globalMatrix;
    //System RHS
    scalarVectorPtr_Type M_globalRightHandSide;

    //System solution (flux and concentration)
    scalarVectorPtr_Type M_fluxAndConcentration;

    // mass matrix
    sparseMatrixPtr_Type M_mediumMassMatrix;

    // mass matrices for fractures
    sparseMatrixPtrContainer_Type M_fractureMassMatrix;

    // System solution (concentration)
    scalarVectorPtr_Type M_mediumConcentrationOld;
    scalarVectorPtr_Type M_mediumConcentration;

    // System fractures solution (velocity, pressure)
    scalarVectorPtrContainer_Type M_fractureConcentrationOld;
    scalarVectorPtrContainer_Type M_fractureConcentration;


    scalarVectorPtr_Type M_stabilizationRightHandSide;

};

typedef TransportFractured TransportFractured_Type;
typedef boost::shared_ptr<TransportFractured_Type> TransportFracturedPtr_Type;

#endif /* TRANSPORTFRACTURED_H_ */
