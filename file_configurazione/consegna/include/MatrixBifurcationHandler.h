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


#ifndef __MATRIXBIFURCATIONHANDLER_H__
#define __MATRIXBIFURCATIONHANDLER_H__ 

#include "Core.h"
#include "Geometry.h"
#include "Parser.h"
#include "FractureHandler.h"
#include <eigen3/Eigen/Dense>

/**************************************************************************/
/*  MatrixBifurcationHandler.h											  */
/*  Classe che gestisce la costruzione delle matrici per il triangolo di  */
/*  intersezione                 										  */
/**************************************************************************/


class MatrixBifurcationHandler
{
public:
	
	/**
	 * Costruttore vuoto di defult. 
	 * Definisce la matrice 3x3 M_Pc costante pari a 1/3 e riempie la matrice 2x2 M_K, leggendo dal file "data".
	 * \param GetPot& dataFile: nome del file data da cui leggere i dati
	 * \param std::string& section = "mediumData/": nome della sezione nel file data in cui leggere, se non è fornita è posta di default 
	 *    											pari a "mediumData/"
	 * \param std::string& subsection = "darcy/": nome della sottosezione nel file data in cui leggere, se non è fornita è posta di default 
	 *    											pari a "darcy/"
	 */ 
	MatrixBifurcationHandler( const GetPot& dataFile,
							  const std::string& type = "Bifurcation",
							  const std::string& section = "mediumData/",
							  const std::string& subsection = "darcy/" );
	
	
	/**
	 * Costruttore che prende in input la matrice di permeabilità.
	 * Definisce la matrice 3x3 M_Pc costante pari a 1/3 e riempie la matrice 2x2 M_K usando la matrice K fornita.
	 * \param Matrix2d K: matrice di permeabilità da usare per modificare M_K
	 */ 
	MatrixBifurcationHandler( Matrix2d K, const std::string& type = "Bifurcation" );
	
	
	/**
	 * Costruttore a partire dalla permeabilità. 
	 * Definisce la matrice 3x3 M_Pc costante pari a 1/3 e riempie la matrice 2x2 M_K usando lo scalare K, inizializzandola come una matrice 
	 * con K sulla diagonale e 0 fuori.
	 */ 
	explicit MatrixBifurcationHandler( scalar_type K, const std::string& type = "Bifurcation" );
	
	
	
	/**
	 * Funzione che inizializza M_intersection, la classe che rappresenta il triangolo di intersezione, e calcola la matrice M_T a partire 
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
	void setK( Matrix2d K ){ M_K=K; }
	
	
	/**
	 * Funzione che restituisce la matrice di permeabilità.
	 * \return Matrix2d M_K: matrice di permeabilità
	 */
	Matrix2d K()const{ return M_K;}
	
	
	/**
	 * Funzione che calcola la matrice M_N, matrice che ha per righe le normali ai lati.
	 * Ovviamente per poter chiamare questo metodo è necessario che il triangolo che rappresenta l'intersezione sia stato già definito.
	 */
	void computeN();
	
	
	/**
	 * Funzione che restituisce la matrice M_N, matrice che ha per righe le normali ai lati.
	 * \return Matrix32 M_N: matrice 3x2
	 */
	Matrix42 N()const { return M_N; }
	

	/**
	 * Funzione che calcola la matrice M_C, matrice che ha per righe i vettori che uniscono il baricentro del triangolo con il punto medio 
	 * di ogni lato.
	 * Ovviamente per poter chiamare questo metodo è necessario che il triangolo che rappresenta l'intersezione sia stato già definito.
	 */
	void computeC();
	
	
	/**
	 * Funzione che restituisce la matrice M_C, matrice che ha per righe i vettori che uniscono il baricentro del triangolo con il punto medio 
	 * di ogni lato.
	 * \return Matrix32 M_C: matrice 3x2
	 */

	Matrix42 C(){ return M_C; }
	
	
	/**
	 * Funzione che costruisce la matrice M_Qc, matrice che ha per colonne una base ortonormale per la matrice trasposta di C.
	 * La funzione richiama il calcolo di M_C, quindi non è necessario che sia stata precedentemente calcolata.
	 */
	void computeQc();
	
	
	/**
	 * Funzione che restituisce la matrice M_Qc, matrice che ha per colonne una base ortonormale per la matrice trasposta di C.
	 * \return Matrix32 M_Qc: matrice 3x2
	 */
	Matrix42 Qc()const { return M_Qc; }
	
		
	/**
	 * Funzione che restituisce la matrice M_Pc, matrice di proiezione sullo spazio nullo di M_Qc.
	 * \return Matrix3d M_Pc: matrice 3x3
	 */
	Matrix4d Pc()const { return M_Pc; }
	
	
	/**
	 * Funzione che calcola la matrice M_T, matrice di trasmissibilità. La matrice M_T deriva dall'applicazione di uno schema alle 
	 * differenze finite mimetiche sul triangolo di intersezione per la legge di Darcy, ignorando l'effetto della gravità.
	 * Questa funzione richiama il calcolo delle matrici M_C, M_N e M_Qc.
	 * L'approssimazione di M_T assume la seguente forma:
	 * 					 T= N K N^T + t \operatorname{trace}(NKN^T)P_c
	 */
	void computeT( scalar_type t = 6.0 );
	
	
	/**
	 * Funzione che calcola la matrice M_T, matrice di trasmissibilità. La matrice M_T deriva dall'applicazione di uno schema alle 
	 * differenze finite mimetiche sul triangolo di intersezione per la legge di Darcy, ignorando l'effetto della gravità.
	 * Questa funzione richiama il calcolo delle matrici M_C, M_N e M_Qc.
	 * L'approssimazione di M_T assume la seguente forma:
	 * 					 T= N K N^T + t \operatorname{trace}(NKN^T)P_c
	 */	
	void computeTsimple( scalar_type t=6.0 );

	
	void computeScap( scalar_type& s, scalar_type t=6.0 );
	
	
	/**
	 * Funzione che inverte una matrice 2x2.
	 * \param Matrix2d& invk: matrice 2x2 da invertire
	 * \param scalar_type det: determinante della matrice
	 */
	void inversion2X2 ( const Matrix2d& invk, scalar_type det );

	/**
	* 
	* 
	* 
	*/	
	void SetDOFIntersecton( FractureHandlerPtr_Type& fractures, scalar_type& DOF );
	
	
	// Operatore di assegnamento
	MatrixBifurcationHandler & operator =( const MatrixBifurcationHandler & mat )
	{
		M_intersection = mat.M_intersection;
		M_K = mat.M_K;
		M_N = mat.M_N;
		M_C = mat.M_C;
		M_Qc = mat.M_Qc;
		M_Pc = mat.M_Pc;
		M_T = mat.M_T;
		
		return *this;
	}

	
	/**
	 * Funzione che restituisce la matrice M_T, matrice di trasmissibilità.
	 * \return Matrix3d M_T, matrice 3x3
	 */
	Matrix4d T()const { return M_T; }
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

private:

	Matrix2d M_K;
	Intersection_Type M_intersection;
	Matrix42 M_N;
	Matrix42 M_C;
	Matrix42 M_Qc;
	Matrix4d M_Pc;
	Matrix4d M_T;
	
	std::string M_type;
};

typedef MatrixBifurcationHandler MatrixBifurcationHandler_Type;									/*!< Classe  MatrixBifurcationHandler */
typedef boost::shared_ptr<MatrixBifurcationHandler_Type> MatrixBifurcationHandlerPtr_Type;		/*!< Puntatore alla classe  MatrixBifurcationHandler */


#endif /* MATRIXBIFURCATIONHANDLER_H_ */