/*
 * TransportStabilization.h
 *
 *  Created on: May 17, 2011
 *      Author: fumagalli
 */

#ifndef TRANSPORTSTABILIZATION_H_
#define TRANSPORTSTABILIZATION_H_ 1

#include "./Core.h"
#include "Parser.h"
#include "MeshHandler.h"
#include "XFEMOperators.h"
#include <algorithm>

class TransportStabilization
{

public:

    TransportStabilization ( const GetPot& dataFile,
                             const std::string& section = "stabilization/" );

    // Stabilization function in a given element
    scalar_type function ( const size_type& iElem ) const;

    void computeMediumPeclet ( const scalarVectorPtr_Type& advectionAbscissa,
                               const scalarVectorPtr_Type& advectionOrdinate,
                               const scalarVector_Type& invDiffusion,
                               const MeshHandlerPtr_Type& mesh,
                               const sparseMatrixPtr_Type& globalMatrix );

    void computeEdgePeclet ( const scalarVector_Type& advection,
                             const scalarVector_Type& invDiffusion,
                             const MeshHandlerPtr_Type& mesh );

    inline const scalarVectorPtr_Type& getMediumPeclet ( ) const
    {
        return M_mediumPeclet;
    }

    inline const scalar_type& getEdgePeclet ( const size_type& dof ) const
    {
        return (*M_edgePeclet) [ dof ];
    }

private:

    // Data for stabilization

    //! Section in the data file
    std::string M_section;

    //! Relazation parameter for the Peclet number
    scalar_type M_relaxation;

    //! Peclet number for each medium elements, without the relaxation parameter
    scalarVectorPtr_Type M_mediumPeclet;

    //! Peclet number for each edges, without the relaxation parameter
    scalarVectorPtr_Type M_edgePeclet;

    //! Stabilization function
    std::string M_stabilizationFunction;

    //! Parser for function evaluation
    mutable LifeV::Parser M_parser;

};

typedef TransportStabilization TransportStabilization_Type;
typedef boost::shared_ptr<TransportStabilization_Type>
        TransportStabilizationPtr_Type;

#endif /* TRANSPORTSTABILIZATION_H_ */
