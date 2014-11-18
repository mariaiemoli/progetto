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


#ifndef EXPORTER_H_
#define EXPORTER_H_ 1

#include "Core.h"
#include "MeshHandler.h"
#include <fstream>


/**************************************************************************/
/*  Exporter.h															  */
/*  Classe in cui definisco le funzioni per esportare i dati              */
/**************************************************************************/


class Exporter
{
public:
	

	/** 
	 * Costruttore che fissa la cartella dove esportare i risultati.
	 * \param dataFile: variabile di tipo GetPot che contiene il nome del file data 
	 * \param section: stringa con il nome della sezione in cui leggere
	 */
    Exporter ( const GetPot& dataFile, const std::string& section = "" );

    
    /**
     * \return std::string& M_vtkFolder: stringa col nome della cartella di esportazione
     */
    inline const std::string& getFolder ( ) const
    {
        return M_vtkFolder;
    }

    
    /**
     * Funzione che esporta una matrice di tipo sparso.
     * \param sparseMatrixPtr_Type& matrix: matrice sparsa
     * \param td::string& nameFile: nome del file in cui esportare 
     */
    void spy ( const sparseMatrixPtr_Type& matrix, const std::string& nameFile ) const;

    
    /**
     * Funzione che esporta un vettore di scalari.
     * \param scalarVectorPtr_Type& vector: vettore di scalari
     * \param std::string& nameFile: nome del file in cui esportare 
     */
    void spy ( const scalarVectorPtr_Type& vector, const std::string& nameFile ) const;

    
    /**
     * Funzione che esporta le regioni della mesh.
     * \param getfem::mesh& mesh: mesh di supporto
     * \param std::string& nameFile = "RegionMesh.vtk": nome del file su cui scrivere la soluzione, se non è dato è posto di default pari a "RegionMesh.vtk"
     */
    void meshRegion ( const getfem::mesh& mesh,
                      const std::string& nameFile = "RegionMesh.vtk" ) const;

private:

    const std::string M_vtkFolder;

};

typedef Exporter Exporter_Type;									/*!< classe Exporter */
typedef boost::shared_ptr<Exporter_Type> ExporterPtr_Type;		/*!< puntatore alla classe Exporter */

#endif /* EXPORTER_H_ */
