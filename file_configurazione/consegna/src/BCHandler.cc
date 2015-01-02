
#include "../include/BCHandler.h"

/**************************************************************************/
/*  BCHandler.cc														  */
/*  Libreria che gestisce le condizioni al bordo nella mesh               */
/**************************************************************************/

BCHandler::BCHandler ( const BCPtrContainer_Type& fractureBC ) :
                       M_fractureBC(fractureBC)
{
}