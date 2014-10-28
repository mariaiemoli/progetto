/**
 * MatrixBifurcationHandler.h
 * 
 * Classe che gestisce la costruzione delle matrici per il triangolo di intersezione
 * 
 */

#ifndef __MATRIXBIFURCATIONHANDLER_H__
#define __MATRIXBIFURCATIONHANDLER_H__1

#include "Core.h"
#include "TriangleHandler.h"
#include "Parser.h"
#include "FractureHandler.h"
#include <eigen3/Eigen/Dense>


class MatrixBifurcationHandler
{
public:
	
	/**
	 * Costruttore vuoto di defult. 
	 * Definisce la matrice 3x3 Pc_ costante pari a 1/3 e riempie la matrice 2x2 K_, leggendo dal file "data"
	 * \param GetPot& dataFile: nome del file data da cui leggere i dati
	 * \param std::string& section = "mediumData/": nome della sezione nel file data in cui leggere, se non è fornita è posta di default 
	 *    											pari a "mediumData/"
	 * \param std::string& subsection = "darcy/": nome della sottosezione nel file data in cui leggere, se non è fornita è posta di default 
	 *    											pari a "darcy/"
	 */ 
	MatrixBifurcationHandler( const GetPot& dataFile,
							  const std::string& section = "mediumData/",
							  const std::string& subsection = "darcy/");
	
	
	/**
	 * Costruttore che prende in input la matrice di permeabilità.
	 * Definisce la matrice 3x3 Pc_ costante pari a 1/3 e riempie la matrice 2x2 K_ usando la matrice K fornita.
	 * \param Matrix2d K: matrice di permeabilità da usare per modificare K_
	 */ 
	MatrixBifurcationHandler( Matrix2d K );
	
	
	/**
	 * Costruttore a partire dalla permeabilità. 
	 * Definisce la matrice 3x3 Pc_ costante pari a 1/3 e riempie la matrice 2x2 K_ usando lo scalare K, inizializzandola come una matrice 
	 * con K sulla diagonale e 0 fuori.
	 */ 
	explicit MatrixBifurcationHandler( scalar_type K );
	
	
	
	/**
	 * Funzione che inizializza M_intersection, la classe che rappresenta il triangolo di intersezione, e calcola la matrice T_ a partire 
	 * dall'insieme delle fratture che costituiscono l'intersezione.
	 * \param FracturePtrContainer_Type& fractures: vettore delle fratture che costituiscono l'intersezione
	 */
	void setMatrices ( FracturePtrContainer_Type& fractures );
	
	
	/**
	 * Funzione che costruisce il triangolo di intersezione a partire da un trangolo già esistente.
	 * \param TriangleData const & T: triangolo già esistente da cui partire per costruire l'intersezione
	 */
	void setTriangle( TriangleData const & T )
	{
		M_intersection.setTriangle ( T );
	}
	
	
	/**
	 * Funzione che cambia la matrice di permeabilità a partire da una matrice data.
	 * \param Matrix2d K: matrice di permeabilità da sostituire
	 */
	void setK( Matrix2d K ){ K_=K; }
	
	
	/**
	 * Funzione che restituisce la matrice di permeabilità.
	 * \return Matrix2d K_: matrice di permeabilità
	 */
	Matrix2d K()const{ return K_;}
	
	
	/**
	 * Funzione che calcola la matrice N_, matrice che ha per righe le normali ai lati.
	 * Ovviamente per poter chiamare questo metodo è necessario che il triangolo che rappresenta l'intersezione sia stato già definito.
	 */
	void computeN();
	
	
	/**
	 * Funzione che restituisce la matrice N_, matrice che ha per righe le normali ai lati.
	 * \return Matrix32 N_: matrice 3x2
	 */
	Matrix32 N()const { return N_; }
	

	/**
	 * Funzione che calcola la matrice C_, matrice che ha per righe i vettori che uniscono il baricentro del triangolo con il punto medio 
	 * di ogni lato.
	 * Ovviamente per poter chiamare questo metodo è necessario che il triangolo che rappresenta l'intersezione sia stato già definito.
	 */
	void computeC();
	
	
	/**
	 * Funzione che restituisce la matrice C_, matrice che ha per righe i vettori che uniscono il baricentro del triangolo con il punto medio 
	 * di ogni lato.
	 * \return Matrix32 C_: matrice 3x2
	 */

	Matrix32 C(){ return C_; }
	
	
	/**
	 * Funzione che costruisce la matrice Qc_, matrice che ha per colonne una base ortonormale per la matrice trasposta di C.
	 * La funzione richiama il calcolo di C_, quindi non è necessario che sia stata precedentemente calcolata.
	 */
	void computeQc();
	
	
	/**
	 * Funzione che restituisce la matrice Qc_, matrice che ha per colonne una base ortonormale per la matrice trasposta di C.
	 * \return Matrix32 Qc_: matrice 3x2
	 */
	Matrix32 Qc()const { return Qc_; }
	
		
	/**
	 * Funzione che restituisce la matrice Pc_, matrice di proiezione sullo spazio nullo di Qc_.
	 * \return Matrix3d Pc_: matrice 3x3
	 */
	Matrix3d Pc()const { return Pc_; }
	
	
	/**
	 * Funzione che calcola la matrice T_, matrice di trasmissibilità. La matrice T_ deriva dall'applicazione di uno schema alle 
	 * differenze finite mimetiche sul triangolo di intersezione per la legge di Darcy, ignorando l'effetto della gravità.
	 * Questa funzione richiama il calcolo delle matrici C_, N_ e Qc_.
	 * L'approssimazione di T_ assume la seguente forma:
	 * 					 T= N K N^T + t \operatorname{trace}(NKN^T)P_c
	 */
	void computeT( scalar_type t = 6.0 );
	
	
	/**
	 * Funzione che calcola la matrice T_, matrice di trasmissibilità. La matrice T_ deriva dall'applicazione di uno schema alle 
	 * differenze finite mimetiche sul triangolo di intersezione per la legge di Darcy, ignorando l'effetto della gravità.
	 * Questa funzione richiama il calcolo delle matrici C_, N_ e Qc_.
	 * L'approssimazione di T_ assume la seguente forma:
	 * 					 T= N K N^T + t \operatorname{trace}(NKN^T)P_c
	 */	
	void computeTsimple( scalar_type t=6.0 );

	
	/**
	 * Funzione che inverte una matrice 2x2.
	 * \param Matrix2d& invk: matrice 2x2 da invertire
	 * \param scalar_type det: determinante della matrice
	 */
	void inversion2X2 ( const Matrix2d& invk, scalar_type det );
	
	
	// Operatore di assegnamento
	MatrixBifurcationHandler & operator =( const MatrixBifurcationHandler & mat )
	{
		M_intersection = mat.M_intersection;
		K_ = mat.K_;
		N_ = mat.N_;
		C_ = mat.C_;
		Qc_ = mat.Qc_;
		Pc_ = mat.Pc_;
		T_ = mat.T_;
	}

	
	/**
	 * Funzione che restituisce la matrice T_, matrice di trasmissibilità.
	 * \return Matrix3d T_, matrice 3x3
	 */
	Matrix3d T()const { return T_; }
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

private:

	Matrix2d K_;
	Intersection_Type M_intersection;
	Matrix32 N_;
	Matrix32 C_;
	Matrix32 Qc_;
	Matrix3d Pc_;
	Matrix3d T_;
};

typedef MatrixBifurcationHandler MatrixBifurcationHandler_Type;									/*!< Classe  MatrixBifurcationHandler */
typedef boost::shared_ptr<MatrixBifurcationHandler_Type> MatrixBifurcationHandlerPtr_Type;		/*!< Puntatore alla classe  MatrixBifurcationHandler */


#endif /* MATRIXBIFURCATIONHANDLER_H_ */