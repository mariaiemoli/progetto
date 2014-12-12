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


#ifndef FRACTUREHANDLER_H_
#define FRACTUREHANDLER_H_ 1

#include "StringUtility.h"
#include "FractureData.h"
#include "BCHandler.h"
#include "LevelSetHandler.h"

/**************************************************************************/
/*  FractureHandler.h													  */
/*  Classe che inizializza e gestisce una frattura             			  */
/**************************************************************************/


class FractureHandler
{
public:

	/**
	 * Flags che uso per distinguere sulla mesh di ogni level set le regioni " non tagliate ", 
	 * cioè dove non vi è intersezione con altre fratture, dalle regioni " tagliate ", cioè dove il level set ha un'intersezione.
	 */
    enum
    {
        FRACTURE_UNCUT = 10000,
        FRACTURE_INTERSECT = 10000
    };

    
    FractureHandler ( const GetPot& dataFile,
                      const size_type& ID,
                      const std::string& section = "fractureData/" );

    /**
     * Funzione che inizializza la frattura.
     * Ogni frattura ha una rappresentazione 1d, ma la mesh reale è una mesh di punti 2d ( sono i punti del tipo (x,y) che stanno sulla 
     * curva y=f(x) ) e questo comporta dei problemi per l'integrazione. Per risolvere il problema per ogni frattura costruiamo una mesh
     * 1d ottenuta proiettando lungo l'asse delle ascisse la mesh reale. Una volta risolto il problema si ritorna sulla mesh reale 
     * interpolando i risultati ottenuti.
     */
    void init ( );


    /**
     * Funzione che calcola il vettore normale alla frattura e la mappa di conversione dalla mesh reale alla mesh piatta.
     */
    void normalVectorAndMap ( const getfem::mesh_fem& mediumMeshFEMPressure );


    /**
     * Funzione che calcola il passo di griglia, h^-1.
     */
    void computeInvH ( const BCHandlerPtr_Type& bcHandler );


    inline const getfem::mesh& getMeshFlat ( ) const
    {
        return M_meshFlat;
    }


    inline getfem::mesh& getMeshFlat ( )
    {
        return M_meshFlat;
    }


    inline const getfem::mesh_fem& getMeshFEMVelocity ( ) const
    {
        return M_meshFEMVelocity;
    }


    inline const getfem::mesh_fem& getMeshFEMPressure ( ) const
    {
        return M_meshFEMPressure;
    }


    inline const getfem::mesh_im& getIntegrationMethodVelocity ( ) const
    {
        return M_integrationMethodVelocity;
    }


    inline const getfem::mesh_im& getIntegrationMethodPressure ( ) const
    {
        return M_integrationMethodPressure;
    }


    inline FractureData& getData ( )
    {
        return M_data;
    }


    inline LevelSetHandlerPtr_Type& getLevelSet ( )
    {
        return M_levelSet;
    }

    inline sizeVector_Type& getDofIntersection( )
    {
        return  M_DOF_Intersection;
    }
	
    inline sizeVector_Type& getDOFBifurcation( )
    {
        return  M_DOF_Bifurcation;
    }
	
    inline void clearDofIntersection( )
    {
        M_DOF_Intersection.clear();
		
		return;
    }
	
    void pushDOFIntersection( size_type i )
    {
        M_DOF_Intersection.push_back( i );
		
		return;
    }

    inline const scalarVector_Type& getEtaNormalInterpolated ( ) const
    {
        return M_etaNormalInterpolated;
    }


    inline const scalarVector_Type& getEtaTangentialInterpolated ( ) const
    {
        return M_etaTangentialInterpolated;
    }


    inline const scalar_type& getEtaTangentialInterpolated ( const size_type& dof ) const
    {
        return M_etaTangentialInterpolated [ dof ];
    }


    inline const scalarVector_Type& getInverseMeshSize ( ) const
    {
        return M_inverseMeshSize;
    }


    inline const scalar_type& getInverseMeshSize ( const size_type& dof ) const
    {
        return M_inverseMeshSize [ dof ];
    }


    inline const scalarVector_Type& getMagnificationMapFactor1 ( ) const
    {
        return M_magnificationMapFactor1;
    }


    inline const scalar_type& getMagnificationMapFactor1 ( const size_type& dof ) const
    {
        return M_magnificationMapFactor1 [ dof ];
    }


    inline const scalarVector_Type& getMagnificationMapFactor2 ( ) const
    {
        return M_magnificationMapFactor2;
    }


    inline const scalar_type& getMagnificationMapFactor2 ( const size_type& dof ) const
    {
        return M_magnificationMapFactor2 [ dof ];
    }


    inline const bgeot::pgeometric_trans& getGeometricTransformation ( ) const
    {
        return M_geometricTransformation;
    }


    inline const scalarVector_Type& getNormal1 ( ) const
    {
        return M_normal1;
    }


    inline const scalarVector_Type& getNormal2 ( ) const
    {
        return M_normal2;
    }


    inline const getfem::mesh& getMeshMapped ( ) const
    {
        return M_meshMapped;
    }


    inline getfem::mesh& getMeshMapped ( )
    {
        return M_meshMapped;
    }

    inline const getfem::mesh_im& getIntegrationMethodPressureVisualization ( ) const
    {
        return M_integrationMethodPressureVisualization;
    }


    inline const getfem::mesh_fem& getMeshFEMPressureVisualization ( ) const
    {
        return M_meshFEMPressureVisualization;
    }


    inline const getfem::mesh_fem& getMeshFEMLinear ( ) const
    {
        return M_meshFEMLinear;
    }


    inline const size_type& getId() const
    {
        return M_ID;
    }


    void numFractures ( const size_type& numFractures )
    {
        M_meshLevelSetIntersect.resize ( numFractures );
        M_levelSetIntersect.resize ( numFractures );
        M_fractureIntersectElements.resize ( numFractures );
        M_fractureIntersectElementsGlobalIndex.resize ( numFractures );
    }

    /**
     * Funzione che, data un'altra frattura con cui si interseca, imposta i  valori legati all'intersezione: costruisce sulla mesh 
     * la regione " tagliata ", aggiunge i gradi di libertà estesi,  e imposta che tali gradi di libertà siano gli stessi sulle due fratture.
     */
    size_type setMeshLevelSetFracture ( FractureHandler& otherFracture, size_type& globalIndex, const std::string& type );


    /**
     * Funzione che restituisce il numero di gradi di libertà estesi per la pressione per la frattura.
     * \return M_extendedPressure.size(): numero di gradi di libertà estesi per la pressione per la frattura.
     */
    size_type getNumExtendedPressure () const
    {
        return M_extendedPressure.size();
    }


    /**
     * Funzione che restituisce ilvettore dei gradi di libertà estesi per la pressione per la frattura.
     * \return M_extendedPressure: vettore dei gradi di libertà estesi per la pressione per la frattura. 
     */
    const sizeVector_Type& getExtendedPressure () const
    {
        return M_extendedPressure;
    }


    /**
     * Funzione che restituisce il numero di gradi di libertà estesi per la velocità per la frattura.
     * \return M_extendedVelocity.size(): numero di gradi di libertà estesi per la velocità per la frattura.
     */
    size_type getNumExtendedVelocity () const
    {
        return M_extendedVelocity.size();
    }


    /**
     * Funzione che restituisce il vettore dei gradi di libertà estesi per la velocità per la frattura.
     * \return M_extendedVelocity: vettore dei gradi di libertà estesi per la velocità per la frattura. 
     */
    const sizeVector_Type& getExtendedVelocity () const
    {
        return M_extendedVelocity;
    }


    /**
     * \return M_meshLevelSetIntersect[f]: mesh del level set di indice f intersecato dalla frattura corrente 
     */
    GFMeshLevelSetPtr_Type getMeshLevelSetIntersect ( const size_type& f )
    {
        return M_meshLevelSetIntersect[f];
    }


    /**
     * \return M_LevelSetIntersect[f]: level set di indice f intersecato dalla frattura corrente 
     */
    GFLevelSetPtr_Type getLevelSetIntersect ( const size_type& f )
    {
        return M_levelSetIntersect[f];
    }


    /**
     * Funzione che calcola in numero totale di intersezione della frattura corrente
     */
    size_type getNumIntersections () const;


    const sizeVectorContainer_Type& getFractureIntersectElements () const
    {
        return M_fractureIntersectElements;
    }


    const pairSizeVectorContainer_Type& getFractureIntersectElementsGlobalIndex () const
    {
        return M_fractureIntersectElementsGlobalIndex;
    }


    pairSizeVectorContainer_Type& getFractureIntersectElementsGlobalIndex ()
    {
        return M_fractureIntersectElementsGlobalIndex;
    }


private:

    size_type M_ID;

    FractureData M_data;

    LevelSetHandlerPtr_Type M_levelSet;
	
	// Salviamo, in caso di una biforcazione, in quale degli estremi della frattura cade
	sizeVector_Type M_DOF_Intersection;
	
	// Nel caso di una biforcazione con due fratture salva il DOF di quella tagliata in due parti
	sizeVector_Type M_DOF_Bifurcation;

    // M_mediummesh per la fratture: M_meshFlat è " piatta " (1d), M_meshMapped è mappata (x(t),y(t))
    getfem::mesh M_meshFlat;
    getfem::mesh M_meshMapped;

    GFMeshLevelSetPtrContainer_Type M_meshLevelSetIntersect;
    GFLevelSetPtrContainer_Type M_levelSetIntersect;

    // eta_gamma = d/K_normale - vector
    scalarVector_Type M_etaNormalInterpolated;
    // eta_t=1/(K_t*d) - vector
    scalarVector_Type M_etaTangentialInterpolated;

    // integration method (velocity)
    getfem::mesh_im M_integrationMethodVelocity;
    // integration method (pressure)
    getfem::mesh_im M_integrationMethodPressure;
    // integration method (pressure)
    getfem::mesh_im M_integrationMethodPressureVisualization;

    //elementi finiti
    // mesh_fem for pressure
    getfem::mesh_fem M_meshFEMPressure;
    // mesh_fem for pressure
    getfem::mesh_fem M_meshFEMPressureVisualization;
    getfem::mesh_fem M_meshFEMVelocity;

    // Fracture extended
    sizeVector_Type M_extendedPressure;
    sizeVector_Type M_extendedVelocity;

    sizeVectorContainer_Type M_fractureIntersectElements;
    pairSizeVectorContainer_Type M_fractureIntersectElementsGlobalIndex;

    //fattore di scala nella mappatura fra frattura piana e mappata
    scalarVector_Type M_magnificationMapFactor1;
    //fattore di scala nella mappatura fra frattura piana e mappata
    scalarVector_Type M_magnificationMapFactor2;

    //componenti della normale
    scalarVector_Type M_normal1;
    scalarVector_Type M_normal2;

    // M_mediumInverseMeshSize = 1.0 / M_mediumMeshSize; per la frattura
    scalarVector_Type M_inverseMeshSize;

    // Geometric transformation usign pressure finite elements type
    bgeot::pgeometric_trans M_geometricTransformation;

    // integration method (velocity)
    getfem::mesh_im M_integrationMethodLinear;
    // mesh_fem for coefficients
    getfem::mesh_fem M_meshFEMLinear;

};

typedef FractureHandler FractureHandler_Type;											/*!< classe FractureHandler */
typedef boost::shared_ptr<FractureHandler_Type> FractureHandlerPtr_Type;				/*!< puntatore alla classe FractureHandler */
typedef std::vector<FractureHandlerPtr_Type> FracturePtrContainer_Type;					/*!< vettore di puntatori alla classe FractureHandler */
typedef std::vector<FractureHandler_Type> FractureContainer_Type;						/*!< vettore di classi FractureHandler */
typedef boost::shared_ptr<FracturePtrContainer_Type> FracturePtrContainerPtr_Type;		/*!< puntatore a un vettore di puntatori alla classe FractureHandler */

#endif /* FRACTUREHANDLER_H_ */
