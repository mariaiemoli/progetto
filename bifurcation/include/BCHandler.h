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

#ifndef BCHANDLER_H_
#define BCHANDLER_H_ 1

#include "Core.h"
#include "BC.h"
#include "FractureData.h"

/**************************************************************************/
/*  BCHandler.h															  */
/*  Classe che gestisce le condizioni al bordo nella mesh                 */
/**************************************************************************/

class BCHandler
{
public:

    // Costruttore
    BCHandler ( const BCPtrContainer_Type& fractureBC );

     
    inline const BCPtr_Type& getFractureBC ( const size_type& id ) const
    {
        return M_fractureBC [ id ];
    }

private:

    BCPtrContainer_Type M_fractureBC;

};

typedef BCHandler BCHandler_Type;											/*!< Classe BCHandler */
typedef boost::shared_ptr<BCHandler_Type> BCHandlerPtr_Type;				/*!< Puntatore alla classe BCHandler */

#endif /* BCHANDLER_H_ */
