/*
 * TimeData.h
 *
 *  Created on: 11/feb/2011
 *      Author: elle
 */

#ifndef TIMEDATA_H_
#define TIMEDATA_H_ 1

#include "Core.h"

class TimeData
{

public:

    TimeData ( const GetPot& dataFile, const std::string& getPotSection = "dataTime/" ) :
        M_getPotSection ( getPotSection ), M_endTime ( dataFile ( ( M_getPotSection + "endTime" ).data (), 1. ) ),
                M_deltaTime ( dataFile ( ( M_getPotSection + "deltaTime" ).data (), 0.1 ) ), M_plotAtEachTimeStep (
                        dataFile ( ( M_getPotSection + "plotAt" ).data (), 1 ) ), M_timeSteps (
                        static_cast<size_type> ( M_endTime / M_deltaTime ) )
    {
    }

    inline scalar_type getEndTime () const
    {
        return M_endTime;
    }

    inline scalar_type getDeltaTime () const
    {
        return M_deltaTime;
    }

    inline size_type getPlotAtEachTimeStep () const
    {
        return M_plotAtEachTimeStep;
    }

    inline size_type getTimeSteps () const
    {
        return M_timeSteps;
    }

private:

    std::string M_getPotSection;
    // tempo finale della simulazione
    scalar_type M_endTime;
    //delta t
    scalar_type M_deltaTime;
    //ogni quanti M_deltaTime plotto la soluzione?!?
    size_type M_plotAtEachTimeStep;
    //numero di step, calcolato con M_endTime e M_deltaTime
    size_type M_timeSteps;

};

typedef TimeData TimeData_Type;
typedef boost::shared_ptr<TimeData_Type> TimeDataPtr_Type;

#endif /* TIMEDATA_H_ */
