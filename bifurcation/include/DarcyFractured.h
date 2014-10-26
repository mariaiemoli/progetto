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
#include "MatrixBifurcationHandler.h"
#include <sstream>
//#include <fstream>


class DarcyFractured
{

public:

    DarcyFractured ( const MediumDataPtr_Type& medium,
                     const MeshHandlerPtr_Type& mesh,
                     const BCHandlerPtr_Type& bcHandler,
                     const FracturesSetPtr_Type& fractures,
                     const ExporterPtr_Type& exporter );
    
    /** 
     * 
     * Funzione che costruisce M_mediumMesh, definisce gli elementi finiti e i metodi di integrazione e seleziona i contorni
     * 
     */
    void init ( );

    
    /** 
     * 
     * Funzione che assembla tutte le matrici e il termine noto di destra 
     * 
     * Il sistema che risulta ha la seguente struttura (DARCY):
     *
     * [  A11   A12 ] [ V ]  = [ Bv ]
     * [ -A12'  0   ] [ P ]  = [ Bp ]
     *
     * (V,P = velocità, pressione)
     *
     * dove
     *
     * BV = Mvd * Vdirichlet + Bstress
     * BP = Mpd * Vdirichlet
     * 
     */
    void assembly ( );

    
    /**
     * Funzione che risolve il sistema ottenuto assemblando il problema di Darcy per tutte le fratture.
     * Una volta risolto il sistema esporta i risultati ottenuti per la pressione per ogni frattura in formato vtk 
     */
    void solve ( );

    
    inline const scalarVectorPtr_Type& getFractureVelocity ( const size_type& f ) const
    {
        return M_fractureVelocity [ f ];
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

    // Vettore che rappresenta l'inverso della permeabilità del mezzo
    scalarVectorContainer_Type M_mediumEtaInterpolated;

    // eta_gamma = d/K_normale - vettore - sui punti della M_mediumMesh grande
    scalarVectorContainer_Type M_fractureEtaNormalOnMedium;

    // Global matrix, darcy
    sparseMatrixPtr_Type M_globalMatrix;
    // Termine noto di destra del sistema
    scalarVectorPtr_Type M_globalRightHandSide;
    // Soluzione del sistema ( velocità + pressione )
    scalarVectorPtr_Type M_velocityAndPressure;

    // Soluzione del sistema nelle fratture
    scalarVectorPtrContainer_Type M_fractureVelocity;
    scalarVectorPtrContainer_Type M_fracturePressure;

    sparseMatrixPtr_Type M_normMatrix; // Norm matrix (i.e., the (H1,L2) norm on

};

typedef DarcyFractured DarcyFractured_Type;
typedef boost::shared_ptr<DarcyFractured_Type> DarcyFracturedPtr_Type;

#endif