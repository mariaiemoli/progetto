/**   
 * DarcyFractured.h
 * structure for the Darcy fractured problem
 *
 */

#ifndef _DARCYFRACTURED_
#define _DARCYFRACTURED_ 1

#include "Core.h"
#include "FracturesSet.h"
#include "XFEMOperators.h"
#include "BCHandler.h"
#include "MediumData.h"
#include "Exporter.h"
#include <sstream>


class DarcyFractured
{

public:

    DarcyFractured ( const MediumDataPtr_Type& medium,
                     const MeshHandlerPtr_Type& mesh,
                     const BCHandlerPtr_Type& bcHandler,
                     const FracturesSetPtr_Type& fractures,
                     const ExporterPtr_Type& exporter );
    
    /** 
     * void init ()
     * Build the M_mediumMesh, set finite element and integration methods and select the boundaries.
     * 
     */
    void init ( );

    
    /** void assembly ()
     * 
     * Assemble all the matrices and the right hand side for the Darcy and a piece of the matrix for the transport
     * 
     * Our system has the following structure (DARCY):
     *
     * [  A11   A12 ] [ V ]  = [ Bv ]
     * [ -A12'  0   ] [ P ]  = [ Bp ]
     *
     * (V,P = velocity, pressure)
     *
     * where
     *
     * BV = Mvd * Vdirichlet + Bstress
     * BP = Mpd * Vdirichlet
     * 
     */
    void assembly ( );

    
    
    void solve ( );

    
    inline const scalarVectorPtr_Type& getMediumVelocity ( ) const
    {
        return M_mediumVelocity;
    }

    
    inline const scalarVectorPtr_Type& getFractureVelocity ( const size_type& f ) const
    {
        return M_fractureVelocity [ f ];
    }

    
    inline const scalarVectorPtr_Type& getMediumVelocityInlet ( ) const
    {
        return M_mediumVelocityInlet;
    }

    
    inline const scalarVectorPtr_Type& getMediumVelocityOutlet ( ) const
    {
        return M_mediumVelocityOutlet;
    }

    
    inline const scalarVectorPtr_Type& getMediumVelocityInterpolatedAbscissa ( ) const
    {
        return M_mediumVelocityInterpolatedAbscissa;
    }

    
    inline const scalarVectorPtr_Type& getMediumVelocityInterpolatedOrdinate ( ) const
    {
        return M_mediumVelocityInterpolatedOrdinate;
    }

    
private:

    // Attributes

    // Data medium
    MediumDataPtr_Type M_mediumData;

    // Mesh
    MeshHandlerPtr_Type M_mesh;

    // BC Handler
    BCHandlerPtr_Type M_bcHandler;

    // Fractures
    FracturesSetPtr_Type M_fractures;

    // Exporter
    ExporterPtr_Type M_exporter;

    // Inverse permeability of the medium - vector
    scalarVectorContainer_Type M_mediumEtaInterpolated;

    // eta_gamma = d/K_normale - vector - sui punti della M_mediumMesh grande
    scalarVectorContainer_Type M_fractureEtaNormalOnMedium;

    // Global matrix, darcy
    sparseMatrixPtr_Type M_globalMatrix;
    // System RHS
    scalarVectorPtr_Type M_globalRightHandSide;
    // System solution (velocity and pressure)
    scalarVectorPtr_Type M_velocityAndPressure;

    // System solution on medium
    scalarVectorPtr_Type M_mediumVelocity;
    scalarVectorPtr_Type M_mediumPressure;

    scalarVectorPtr_Type M_mediumVelocityInlet;
    scalarVectorPtr_Type M_mediumVelocityOutlet;

    // System solution on fractures
    scalarVectorPtrContainer_Type M_fractureVelocity;
    scalarVectorPtrContainer_Type M_fracturePressure;

    sparseMatrixPtr_Type M_normMatrix; // Norm matrix (i.e., the (H1,L2) norm on
    
    /**
     * the (velocity,pressure) space is given by 
     * ||M_darcyVelocityAndPressure||^2 = M_darcyVelocityAndPressure' * M_normMatrix * M_darcyVelocityAndPressure, 
     * where M_darcyVelocityAndPressure = (M_darcyMediumVelocity, M_darcyMediumPressure) (see below)
     * 
     */ 

    // Interpolation of the velocity in each baricenter
    scalarVectorPtr_Type M_mediumVelocityInterpolatedAbscissa;
    scalarVectorPtr_Type M_mediumVelocityInterpolatedOrdinate;

};

typedef DarcyFractured DarcyFractured_Type;
typedef boost::shared_ptr<DarcyFractured_Type> DarcyFracturedPtr_Type;

#endif