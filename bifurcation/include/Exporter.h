/** Exporter.h
 *
 * Created on: Apr 13, 2011
 *
 * \author Alessio Fumagalli
 * 
 * libreria in cui definisco le funzioni per esportare i dati
 *
 */

#ifndef EXPORTER_H_
#define EXPORTER_H_ 1

#include "Core.h"
#include "MeshHandler.h"
#include <fstream>

class Exporter
{
public:
	

	/** 
	 * Costruttore che fissa la cartella dove esportare i risultati
	 * \param dataFile: variabile di tipo GetPot che contiene il nome del file data 
	 * \param section: stringa con il nome della sezione in cui leggere
	 */
    Exporter ( const GetPot& dataFile, const std::string& section = "" );

    
    /**
     * \return M_vtkFolder: stringa col nome della cartella di esportazione
     */
    inline const std::string& getFolder ( ) const
    {
        return M_vtkFolder;
    }

    
    /**
     * funzione che esporta una matrice di tipo sparso
     * \param matrix: matrice sparsa
     * \param nameFile: nome del file in cui esportare 
     */
    void spy ( const sparseMatrixPtr_Type& matrix, const std::string& nameFile ) const;

    
    /**
     * funzione che esporta un vettore di scalari
     * \param vector: vettore di scalari
     * \param nameFile: nome del file in cui esportare 
     */
    void spy ( const scalarVectorPtr_Type& vector, const std::string& nameFile ) const;

    
    /**
     * funzione che esporta le regioni della mesh
     * \param mesh: mesh di supporto
     * \param nameFile: nome del file su cui scrivere la soluzione
     */
    void meshRegion ( const getfem::mesh& mesh,
                      const std::string& nameFile = "RegionMesh.vtk" ) const;

private:

    const std::string M_vtkFolder;

};

typedef Exporter Exporter_Type;									/*!< classe Exporter */
typedef boost::shared_ptr<Exporter_Type> ExporterPtr_Type;		/*!< puntatore alla classe Exporter */

#endif /* EXPORTER_H_ */
