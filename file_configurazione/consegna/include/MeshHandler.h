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


#ifndef MESHHANDLER_H_
#define MESHHANDLER_H_ 1

#include "Core.h"
#include "UsefulFunctions.h"
#include "FracturesSet.h"
#include "BCHandler.h"
#include <math.h>

/**************************************************************************/
/*  MeshHandler.h														  */
/*  Classe che costruisce e manipola la mesh di supporto                  */
/**************************************************************************/

class MeshHandler
{
public:

	enum
	{
		UNCUT_REGION = 100
	};


	/** 
	 *
	 * costruttore della classe MeshHandler
	 * riempie i vari campi della classe con i dati nel file data
	 *
	 * la classe contiente informazioni circa le dimensione del dominio, il tipo di mesh,
	 * la discretizzazione e i metodi di integrazione
	 * 
	 * \param dataFile: variabile di tipo GetPot che contiene il nome del file data 
	 * \param section: stringa con il nome della sezione in cui leggere
	 *
	 */
	MeshHandler ( const GetPot& dataFile, const std::string& sectionDomain = "");


	/** 
	 *
	 * funzione che costruisce la mesh
	 * a seconda che il campo "meshExternal" nel file data sia posto pari a "none" o altro
	 * costruisce la mesh usando i dati del file o la importa dall'esterno
	 */
    void setUpMesh ( );


	/** 
	 *
	 * funzione che definisce le regioni della mesh
	 * \param fracture: puntatore all'insieme delle fratture
	 */
    void setUpRegions ( const FracturesSetPtr_Type& fracture );


	/** 
	 *
	 * funzione che controlla se un convesso è effettivamente tagliato da una frattura 
	 * \param nodes: insieme dei nodi geometrici che definiscono il convesso
	 * \param fracture: puntatore a una frattura
	 */
    bool compreso ( const bgeot::basic_mesh::ref_mesh_pt_ct nodes, FractureHandlerPtr_Type& fracture );


    /**
     * funzione che definisce gli elementi finiti su cui si lavora
     *
     */
    void setUpFEM ( );

    /**
     * funzione che calcola l'inverso del passo della mesh, h^-1
     */
    void computeMeshMeasures ( );


    /** 
     *
     * funzione che esporta in formato .vtk per paraview gli elementi che sono tagliati dal level set M_levelSet
     * pone pari a 1 se sono tagliati e 0 se non lo sono
     * \param vtkFolder: nome della cartella in cui esportare
     * \param fileName: nome da dare al file da esportare con la soluzione
     */
    void printCuttedElements ( const std::string& vtkFolder = "vtk/", const std::string& fileName = "CuttedElements" ) const;


    /**
     * \return M_spatialDiscretization: discretizzazione spaziale, numero di elementi in cui ogni lato del dominio viene suddiviso
     */
    inline size_type getSpatialDiscretization ( ) const
    {
        return M_spatialDiscretization;
    }

    /**
     * \return M_inclination: grandezza che indica di quanto il dominio reale è inclinato rispetto al dominio di riferimento
     */
    inline scalar_type getInclination ( ) const
    {
        return M_inclination;
    }

    
    /**
     * \return M_lengthAbscissa: grandezza che indica la lunghezza dell'ascissa del dominio reale
     */
    inline scalar_type getLengthAbscissa ( ) const
    {
        return M_lengthAbscissa;
    }


    /**
     * \return M_lengthOrdinate: grandezza che indica la lunghezza dell'ordinata del dominio reale
     */
    inline scalar_type getLengthOrdinate ( ) const
    {
        return M_lengthOrdinate;
    }


    /**
     * \return M_lengthQuota: grandezza che indica la lunghezza della quota del dominio reale
     */
    inline scalar_type getLengthQuota ( ) const
    {
        return M_lengthQuota;
    }


    /**
     * \return M_meshType: std::string M_meshType, grandezza che rappresenta una trasformazione per GetFEM++, definisce l'elemento di riferimento su cui è 
     * descritto il metodo di integrazione
     */
    inline std::string getMeshType ( ) const
    {
        return M_meshType;
    }


    inline const getfem::mesh& getMesh ( ) const
    {
        return M_mesh;
    }


    inline getfem::mesh& getMesh ( )
    {
        return M_mesh;
    }


    inline getfem::mesh_level_set& getMeshLevelSet ()
    {
        return M_meshLevelSet;
    }


    inline bgeot::dim_type getSpaceDimension ( ) const
    {
        return M_spaceDimension;
    }


    inline std::string getIntegrationTypeVelocity ( ) const
    {
        return M_integrationTypeVector;
    }


    inline std::string getIntegrationTypePressure ( ) const
    {
        return M_integrationTypeScalar;
    }


    inline std::string getFEMTypeVelocity ( ) const
    {
        return M_fEMTypeVector;
    }


    inline std::string getFEMTypePressure ( ) const
    {
        return M_fEMTypeScalar;
    }


    inline const getfem::mesh_fem& getMeshFEMCoefficients ( ) const
    {
        return M_meshFEMCoefficients;
    }


    inline const getfem::mesh_fem& getMeshFEMScalar ( ) const
    {
        return M_meshFEMScalar;
    }


    inline const getfem::mesh_fem& getMeshFEMVector ( ) const
    {
        return M_meshFEMVector;
    }


    inline const getfem::mesh_im& getIntegrationMethodVector ( ) const
    {
        return M_integrationMethodVector;
    }


    inline const getfem::mesh_im& getIntegrationMethodScalar ( ) const
    {
        return M_integrationMethodScalar;
    }


    inline const getfem::mesh_region& getRegion ( const size_type& regionFlag ) const
    {
        return M_mesh.region(regionFlag);
    }


    inline getfem::mesh_region& getRegion ( const size_type& regionFlag )
    {
        return M_mesh.region(regionFlag);
    }


    inline const scalarVector_Type& getMeshSize ( ) const
    {
        return M_meshSize;
    }


    inline const scalarVector_Type& getInverseMeshSize ( ) const
    {
        return M_inverseMeshSize;
    }


    inline const scalar_type& getInverseMeshSizeDOF ( const size_type& dof ) const
    {
        return M_inverseMeshSize [ dof ];
    }


    inline const sizeVector_Type& getExtendedDOFScalar ( const size_type& id ) const
    {
        return M_extendedDOFScalar [ id ];
    }


    /**
     * funzione che conta il numero di gradi di libertà estesi per la pressione per tutte le fratture che hanno indice < id
     * \param id: indice che identifica una frattura
     * \return numero dei gradi di libertà 
     */
    size_type getCountExtendedDOFScalar ( const scalar_type& id ) const;

    
    /**
     * funzione che conta il numero di gradi di libertà estesi per la velocità per tutte le fratture che hanno indice < id
     */
    size_type getCountExtendedDOFVector ( const scalar_type& id ) const;

    
    inline const size_type& getExtendedDOFScalar ( const size_type& id,
                                                   const size_type& dof ) const
    {
        return M_extendedDOFScalar [ id ] [ dof ];
    }


    inline const sizeVector_Type& getExtendedDOFVector ( const size_type& id ) const
    {
        return M_extendedDOFVector [ id ];
    }


    inline const size_type& getExtendedDOFVector ( const size_type& id,
                                                   const size_type& dof ) const
    {
        return M_extendedDOFVector [ id ] [ dof ];
    }


    inline size_type getCountExtendedIntersectDOFScalar () const
    {
        return M_extendedIntersectDOFScalar.size();
    }


    const sizeVector_Type& getExtendedIntersectDOFScalar() const
    {
        return M_extendedIntersectDOFScalar;
    }


    inline size_type getCountExtendedIntersectDOFVector () const
    {
        return M_extendedIntersectDOFVector.size();
    }


    const sizeVector_Type& getExtendedIntersectDOFVector() const
    {
        return M_extendedIntersectDOFVector;
    }


    inline const bgeot::pgeometric_trans& getGeometricTransformation ( ) const
    {
        return M_geometricTransformation;
    }


    inline sizeVector_Type& getNonCut ( )
    {
        return M_nonCut;
    }


    inline const scalar_type& getEdgeLength ( const size_type& dof ) const
    {
        return M_edgeLength [ dof ];
    }


    inline const scalar_type& getCircumcentersDistance ( const size_type& dof ) const
    {
        return M_circumcentersDistance [ dof ];
    }


private:
    
    /*
     * questa sistema la cut region creando una lista dei triangoli che non sono veramente tagliati,
     * in pratica quando A1 o A2 sono minori di una certa tolleranza
     */
    void fixCutRegion ( const FractureHandlerPtr_Type& fracture );


    // the M_mediumMesh
    getfem::mesh M_mesh;

    getfem::mesh_level_set M_meshLevelSet;

    std::string M_meshType;

    bgeot::dim_type M_spaceDimension;


    /**elementi finiti
     *
     * mesh fem for velocity
     *
     * mesh_fem for pressure
     *
     * mesh_fem for coefficients
     *
     */

    std::string M_fEMTypeVector;
    getfem::mesh_fem M_meshFEMVector;

    std::string M_fEMTypeScalar;
    getfem::mesh_fem M_meshFEMScalar;

    getfem::mesh_fem M_meshFEMCoefficients;


    /** metodi di integrazione
     *
     * integration method for vector fields
     *
     * integration method for scalar fields
     *
     */

    std::string M_integrationTypeVector;
    getfem::mesh_im M_integrationMethodVector;

    std::string M_integrationTypeScalar;
    getfem::mesh_im M_integrationMethodScalar;


    // Geometric transformation usign pressure finite elements type
    bgeot::pgeometric_trans M_geometricTransformation;

    // M_mediumMeshSize = local M_mediumMesh size
    scalarVector_Type M_meshSize;
    // M_mediumInverseMeshSize = 1.0 / M_mediumMeshSize;
    scalarVector_Type M_inverseMeshSize;

    std::string M_meshExternal;
    std::string M_meshFolder;

    size_type M_spatialDiscretization;
    scalar_type M_inclination;
    scalar_type M_lengthAbscissa;
    scalar_type M_lengthOrdinate;
    scalar_type M_lengthQuota;

    /**gradi di libertà estesi (cioè raddoppiati nei triangoli tagliati)
     *
     * The extended dofs for p
     * Base, one sizeVector for each fracture
     * Intersect, one sizeVector for each interseciton type
     *
     * The extended dofs for u
     * Base, one sizeVector for each fracture
     * Intersect, one sizeVector for each intersection type
     *
     */

    sizeVectorContainer_Type M_extendedDOFScalar;
    sizeVector_Type M_extendedIntersectDOFScalar;

    sizeVectorContainer_Type M_extendedDOFVector;

    sizeVector_Type M_extendedIntersectDOFVector;

    // The circumcenters of all the elements
    scalarVector_Type M_circumcentersAbscissa;
    scalarVector_Type M_circumcentersOrdinate;
    scalarVector_Type M_circumcentersDistance;

    // The midpoints for each edge
    scalarVector_Type M_edgeMidpointAbscissa;
    scalarVector_Type M_edgeMidpointOrdinate;

    // Edge length
    scalarVector_Type M_edgeLength;

    // per sistemare il caso degenere (triangolo non veramente tagliato)
    sizeVector_Type M_nonCut;

};

typedef MeshHandler MeshHandler_Type;									/*!< Classe MeshHandler */
typedef boost::shared_ptr<MeshHandler_Type> MeshHandlerPtr_Type;		/*!< Puntatore alla classe MeshHandler */


#endif /* MESHHANDLER_H_ */
