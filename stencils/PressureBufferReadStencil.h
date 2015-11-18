#ifndef _PRESSURE_BUFFER_READ_STENCIL_H_
#define _PRESSURE_BUFFER_READ_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include <string>

/**
 *
 * Reads the pressure values at the boundary from an array.
 */
class PressureBufferReadStencil : public BoundaryStencil<FlowField> {

private:
    FLOAT * _pressuresLeft;         // Array for reading the pressures at the left boundary
    FLOAT * _pressuresRight;        // Array for reading the pressures at the right boundary
    FLOAT * _pressuresBottom;       // Array for reading the pressures at the bottom boundary
    FLOAT * _pressuresTop;          // Array for reading the pressures at the top boundary
    FLOAT * _pressuresFront;        // Array for reading the pressures at the front boundary
    FLOAT * _pressuresBack;         // Array for reading the pressures at the back boundary

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
        FLOAT * pressuresBottom, FLOAT * pressuresTop);

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
        FLOAT * pressuresFront, FLOAT * pressuresBack);

    //@brief Functions for the 2D problem. Coordinates entered in alphabetical order.
    //@{
    void applyLeftWall   ( FlowField & flowField, int i, int j );
    void applyRightWall  ( FlowField & flowField, int i, int j );
    void applyBottomWall ( FlowField & flowField, int i, int j );
    void applyTopWall    ( FlowField & flowField, int i, int j );
    //@}

    //@brief Functions for the 3D problem. Coordinates entered in alphabetical order.
    //@{
    void applyLeftWall   ( FlowField & flowField, int i, int j, int k );
    void applyRightWall  ( FlowField & flowField, int i, int j, int k );
    void applyBottomWall ( FlowField & flowField, int i, int j, int k );
    void applyTopWall    ( FlowField & flowField, int i, int j, int k );
    void applyFrontWall  ( FlowField & flowField, int i, int j, int k );
    void applyBackWall   ( FlowField & flowField, int i, int j, int k );
    //@}

};

#endif
