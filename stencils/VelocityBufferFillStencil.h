#ifndef _VELOCITY_BUFFER_FILL_STENCIL_H_
#define _VELOCITY_BUFFER_FILL_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include <string>

/**
 *
 * Fills an array with the velocity values at the boundary.
 */
class VelocityBufferFillStencil : public BoundaryStencil<FlowField> {

private:
    FLOAT * _velocitiesLeft;         // Array for saving the velocities at the left boundary
    FLOAT * _velocitiesRight;        // Array for saving the velocities at the right boundary
    FLOAT * _velocitiesBottom;       // Array for saving the velocities at the bottom boundary
    FLOAT * _velocitiesTop;          // Array for saving the velocities at the top boundary
    FLOAT * _velocitiesFront;        // Array for saving the velocities at the front boundary
    FLOAT * _velocitiesBack;         // Array for saving the velocities at the back boundary
    int _lowOffset;

    //@brief Refactored stencil functionality
    //@{
    inline void applyStencil2D(FlowField & flowField, FLOAT * velBuffer, int i, int j, int indVar);
    inline void applyStencil3D(FlowField & flowField, FLOAT * velBuffer, int i, int j, int k, int dimFast, int indSlow, int indFast);
    //@}

public:

    /** Constructor for 2D case
     * @param parameters Parameters of the simulation
     * @param velocitiesLeft Pointer to an array of length 2*(N+2) for saving the velocities at the left boundary
     * @param velocitiesRight Pointer to an array of length 2*(N+2) for saving the velocities at the right boundary
     * @param velocitiesBottom Pointer to an array of length 2*(N+2) for saving the velocities at the bottom boundary
     * @param velocitiesTop Pointer to an array of length 2*(N+2) for saving the velocities at the top boundary
     */
    VelocityBufferFillStencil(const Parameters & parameters,
        FLOAT * velocitiesLeft, FLOAT * velocitiesRight,
        FLOAT * velocitiesBottom, FLOAT * velocitiesTop,
        int lowOffset = 0);

    /** Constructor for 3D case
     * @param parameters Parameters of the simulation
     * @param velocitiesLeft Pointer to an array of length 3*(N+2)^2 for saving the velocities at the left boundary
     * @param velocitiesRight Pointer to an array of length 3*(N+2)^2 for saving the velocities at the right boundary
     * @param velocitiesBottom Pointer to an array of length 3*(N+2)^2 for saving the velocities at the bottom boundary
     * @param velocitiesTop Pointer to an array of length 3*(N+2)^2 for saving the velocities at the top boundary
     * @param velocitiesFront Pointer to an array of length 3*(N+2)^2 for saving the velocities at the front boundary
     * @param velocitiesBack Pointer to an array of length 3*(N+2)^2 for saving the velocities at the back boundary
     */
    VelocityBufferFillStencil(const Parameters & parameters,
        FLOAT * velocitiesLeft, FLOAT * velocitiesRight,
        FLOAT * velocitiesBottom, FLOAT * velocitiesTop,
        FLOAT * velocitiesFront, FLOAT * velocitiesBack,
        int lowOffset = 0);

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
