#ifndef _PRESSURE_BUFFER_READ_STENCIL_H_
#define _PRESSURE_BUFFER_READ_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include "ScalarBufferReadStencil.h"
#include <string>

/**
 *
 * Reads the pressure values at the boundary from an array.
 */
class PressureBufferReadStencil : public ScalarBufferReadStencil<FlowField> {

public:

    /** Constructor for 2D case
     * @param parameters Parameters of the simulation
     * @param pressuresLeft Pointer to an array of length (N+2) for reading the pressures at the left boundary
     * @param pressuresRight Pointer to an array of length (N+2) for reading the pressures at the right boundary
     * @param pressuresBottom Pointer to an array of length (N+2) for reading the pressures at the bottom boundary
     * @param pressuresTop Pointer to an array of length (N+2) for reading the pressures at the top boundary
     */
    PressureBufferReadStencil(const Parameters & parameters,
        FLOAT * pressuresLeft, FLOAT * pressuresRight,
        FLOAT * pressuresBottom, FLOAT * pressuresTop,
        int lowOffset = 0);

    /** Constructor for 3D case
     * @param parameters Parameters of the simulation
     * @param pressuresLeft Pointer to an array of length (N+2)^2 for reading the pressures at the left boundary
     * @param pressuresRight Pointer to an array of length (N+2)^2 for reading the pressures at the right boundary
     * @param pressuresBottom Pointer to an array of length (N+2)^2 for reading the pressures at the bottom boundary
     * @param pressuresTop Pointer to an array of length (N+2)^2 for reading the pressures at the top boundary
     * @param pressuresFront Pointer to an array of length (N+2)^2 for reading the pressures at the front boundary
     * @param pressuresBack Pointer to an array of length (N+2)^2 for reading the pressures at the back boundary
     */
    PressureBufferReadStencil(const Parameters & parameters,
        FLOAT * pressuresLeft, FLOAT * pressuresRight,
        FLOAT * pressuresBottom, FLOAT * pressuresTop,
        FLOAT * pressuresFront, FLOAT * pressuresBack,
        int lowOffset = 0);

protected:
    FLOAT & getScalar( FlowField & flowField, int i, int j );
    FLOAT & getScalar( FlowField & flowField, int i, int j, int k );

};

#endif
