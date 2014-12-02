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


#ifndef _FRACTUREINTERSECT_
#define _FRACTUREINTERSECT_ 1

#include "Core.h"
#include "IntersectData.h"
#include "FractureHandler.h"
#include "UsefulFunctions.h"
#include <assert.h>
#include <map>

/**************************************************************************/
/*  FractureIntersect.h													  */
/*  Classe che contiene tutte le intersezioni divise in base al loro tipo */
/*  ( Cross, Bifurcation, Parallel )                					  */
/**************************************************************************/


class FractureIntersect
{
public:
        enum IntersectionType
        {
                Parallel = 400000,
                Cross = 500000,
                Bifurcation = 600000,
				Bifurcation2 = 700000
        };


        typedef std::map < IntersectionType, IntersectDataContainer_Type > mapIntersection_Type;
        typedef std::pair < size_type, size_type > regionLevelSetPair_Type;

        
        /**
         * Costruttore nullo.
         * L'idea è quella di assegnare dei flags ad ogni sottoregione per poter distinguere il tipo di intersezione.
         * In particolare abbiamo tre tipi di intersezione:
         * 		- Cross: quando due fratture di intersecano in un punto formando una X
         * 		- Bifurcation: quando tre fratture si intersecano in un punto comune formando una Y
         * 		- Parallel: quando due o più fratture passano nello stesso elemento della mesh di supporto ma non si intersecano, 
         * 		    		almeno non in quell'elemento
         */
        FractureIntersect ();

        
        /**
         * Funzione che costruisce la classe delle intersezioni.
         * Per ogni elemento della mesh di supporto verifica quanti level set lo attraversano. Se vi passano due o più fratture
         * lo considero come una possibile intersezione e verifico se effettivamente lo è e di che tipo.
         * \param GetPot& dataFile: dome del file data dove leggere i dati
         * \param getfem::mesh_level_set& meshLevelSet: mesh che tiene conto della presenza dei levelset e che va aggiornata in
         * 												presenza di intersezioni
         * \param FracturePtrContainer_Type& fractures: puntatore all'insieme di tutte le fratture
         * 
         */
        void constructIntesection ( const getfem::mesh& mesh, getfem::mesh_level_set& meshLevelSet,
                                    const FracturePtrContainer_Type& fractures );

        
        /**
         * Funzione che restituisce tutte le intersezioni del tipo richiesto.
         * \param IntersectionType type: contiene il tipo delle intersezioni che mi interessano
         * \return IntersectDataContainer_Type M_intersections [ type ]: restituisce il vettore di tutte le intersezioni del tipo richiesto
         */
        IntersectDataContainer_Type& getIntersectionsOfType ( IntersectionType type )
        {
                return M_intersections [ type ];
        }
        
        
        /**
         * Funzione che restituisce tutte le intersezioni di tipo " Cross ".
         * \return IntersectDataContainer_Type: restituisce il vettore di tutte le intersezioni del tipo " Cross "
         */
        IntersectDataContainer_Type getCrossIntersections () const;
        

        /**
         * Funzione che restituisce tutte le intersezioni di tipo " Bifurcation ".
         * \return IntersectDataContainer_Type: restituisce il vettore di tutte le intersezioni del tipo " Bifurcation "
         */
        IntersectDataContainer_Type getBifurcationIntersections () const;
		
       /**
        * Funzione che restituisce tutte le intersezioni di tipo " Bifurcation2 ".
        * \return IntersectDataContainer_Type: restituisce il vettore di tutte le intersezioni del tipo " Bifurcation2 "
        */
       IntersectDataContainer_Type getBifurcation2Intersections () const;
        

        /**
         * Funzione che restituisce tutto l'insieme delle intersezioni con i rispettivi tipi ( una mappa ).
         * \return mapIntersection_Type: std::map < IntersectionType, IntersectDataContainer_Type >
         */
        mapIntersection_Type& getIntersections ()
        {
                return M_intersections;
        }


        /** 
         * Funzione  restituisce il numero di intersezioni trovate del tipo type.
         * \param IntersectionType: enum IntersectionType { Parallel = 400000, Cross = 500000, Bifurcation = 600000 }, 
         * 							contiene il tipo delle intersezioni che mi interessano
         */
        size_type getNumberIntersectionOfType ( IntersectionType type ) const;
        
        
        /**
         * Funzione che restituisce il numero di intersezioni trovate di tipo " Cross ".
         */
        size_type getNumberCross () const;
        
        
        /**
		 * Funzione che restituisce il numero di intersezioni trovate di tipo " Bifurcation ".
		 */
        size_type getNumberBifurcation () const;
		
       /**
	    * Funzione che restituisce il numero di intersezioni trovate di tipo " Bifurcation2 ".
	    */
       size_type getNumberBifurcation2 () const;


        /** 
         * Funzione che restituisce il numero totale di intersezioni
         */
        size_type getNumberIntersections () const;


        size_type getNumberType () const;


        size_type getBasisFunctionOfType ( IntersectionType type ) const;
        
		
		void setDOFIntersection( const getfem::mesh& M_mesh, FractureHandlerPtr_Type& fracture , size_type i );
		

		size_type FindDOF_Intersection( const bgeot::basic_mesh::ref_mesh_pt_ct nodes, FractureHandlerPtr_Type& fracture );
		
		bool checkCross( const sizeVector_Type& levelSets, const FracturePtrContainer_Type& fractures );
		
		void isRealIntersection ( const getfem::mesh& M_mesh, const size_type i, sizeVector_Type& levelSet, FracturePtrContainer_Type& fracturesInvolved );



private:


        /**
         * Funzione che restituisce il tipo di intersezione, " Cross ", " Bifurcation " o " Parallel".
         * L'idea è quella di considerare gli elementi della mesh di supporto in cui passano due o più levelset e:
         * - se non vi è intersezione la funzione restituisce " Parallel "
         * - se due levelset si intersecano la funzione restituisce " Cross "
         * - se tre levelset si interesecano la funzione restituisce " Bifurcation "
         * \param getfem::mesh_level_set& meshLevelSet: mesh di supporto con i levelset
         * \param size_type& elementID: ID dell'elemento di supporto dove passano i level set
         * \param sizeVector_Type& levelSets: level set interessati
         */
        IntersectionType intersectionType ( getfem::mesh_level_set& meshLevelSet,
                                            const size_type& elementID,
                                            const sizeVector_Type& levelSets,
											const FracturePtrContainer_Type& fractures );


        scalar_type integrateWithBooleanOperation ( getfem::mesh_level_set& meshLevelSet,
                                                    const size_type& elementID,
                                                    const std::string& operation ) const;


        mapIntersection_Type M_intersections;
        std::map < regionLevelSetPair_Type, IntersectionType > M_subRegionIntersection;
		std::map < IntersectionType, size_type> M_basisFunctionOfType;

};

typedef FractureIntersect FractureIntersect_Type;									/*!< classe FractureIntersect */
typedef boost::shared_ptr < FractureIntersect > FractureIntersectPtr_Type;			/*!< puntatore alla classe FractureIntersect */ 

#endif /* _FRACTUREINTERSECT_H_ */
