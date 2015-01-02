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


#ifndef _INTERSECTDATA_
#define _INTERSECTDATA_ 1

#include "FractureHandler.h"

/**************************************************************************/
/*  IntersectData.h														  */
/*  Classe che contiene tutte le informazioni su un'intersezione		  */
/**************************************************************************/


class IntersectData
{

public:


    IntersectData ()
    {
    } // costruttore nullo
	

    IntersectData ( const IntersectData& in )
    {
    	copy ( in );
    } // costruttore di copia


    
    /**
     * Funzione che assegna a M_fractures le fratture che si intersecano, M_elementID l'ID dell'elemento in cui avviene l'intersezione e,
     * se l'intersezione è di tipo " Bifurcation ", costruisce il triangolo dell'intersezione.
     */
    void setIntersection ( const size_type& elementID,
                           const FracturePtrContainer_Type& fractures )
    {
        M_elementID = elementID;
        M_fractures = fractures;
        
        return;
        
    }

    
    void setDOFPosition ( const sizeVector_Type& dofPressure,
                          const sizeVectorContainer_Type& dofVelocity )
    {
        M_dofPressure = dofPressure;
        M_dofVelocity = dofVelocity;
    } // setDofPosition


    size_type getDOFPressure ( const size_type& component ) const
    {
        return M_dofPressure [ component ];
    } // getDOFPressure


    size_type getDOFVelocity ( const size_type& component, const size_type& region ) const
    {
        return M_dofVelocity [ component ] [ region ];
    } //  getDOFVelocity


    const sizeVector_Type& getDOFVelocity ( const size_type& region ) const
    {
        return M_dofVelocity [ region ];
    } // getDOFVelocity


    const FracturePtrContainer_Type& getFractures () const
    {
        return M_fractures;
    } // getFractures


    const FractureHandlerPtr_Type& getFracture ( const size_type& f ) const
    {
        return M_fractures [ f ];
    } // getFracture


    size_type getNumFractures () const
    {
        return M_fractures.size();
    } // getNumFractures


    const size_type& getElementID () const
    {
        return M_elementID;
    } // getElementID


    const IntersectData& operator = ( const IntersectData & in )
    {
        copy ( in );
        
        return *this;
    } // operator =


private:

    void copy ( const IntersectData& in );

    
    // Insieme delle fratture che si intersecano
    FracturePtrContainer_Type M_fractures;
    
    // ID dell'elemento nella mesh di supporto dove le fratture si intersecani
    size_type M_elementID;
    
    sizeVector_Type M_dofPressure;
    sizeVectorContainer_Type M_dofVelocity;

};

typedef IntersectData IntersectData_Type;											/*!< Classe IntersectData */
typedef std::vector < IntersectData_Type > IntersectDataContainer_Type;				/*!< Vettore di classi IntersectData */
typedef boost::shared_ptr < IntersectData_Type > IntersectDataPtr_Type;				/*!< Puntatore alla classe IntersectData */
typedef std::vector < IntersectDataPtr_Type > IntersectDataPtrContainer_Type;		/*!< Vettore di puntatori alla classe IntersectData */


#endif /* _INTERSECTDATA_H_ */