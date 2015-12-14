#ifndef _SCALAR_BUFFER_FILL_STENCIL_H_
#define _SCALAR_BUFFER_FILL_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include <string>

/**
 *
 * Fills an array with the values at the boundary.
 */
template<class FlowField>
class ScalarBufferFillStencil : public BoundaryStencil<FlowField> {

protected:
    FLOAT * _bufferLeft;         // Array for saving the values at the left boundary
    FLOAT * _bufferRight;        // Array for saving the values at the right boundary
    FLOAT * _bufferBottom;       // Array for saving the values at the bottom boundary
    FLOAT * _bufferTop;          // Array for saving the values at the top boundary
    FLOAT * _bufferFront;        // Array for saving the values at the front boundary
    FLOAT * _bufferBack;         // Array for saving the values at the back boundary
    int _lowOffset;

public:

    /** Constructor for 2D case
     * @param parameters Parameters of the simulation
     * @param bufferLeft Pointer to an array of length (N+2) for saving the values at the left boundary
     * @param bufferRight Pointer to an array of length (N+2) for saving the values at the right boundary
     * @param bufferBottom Pointer to an array of length (N+2) for saving the values at the bottom boundary
     * @param bufferTop Pointer to an array of length (N+2) for saving the values at the top boundary
     */
    ScalarBufferFillStencil(const Parameters & parameters,
            FLOAT * bufferLeft, FLOAT * bufferRight,
            FLOAT * bufferBottom, FLOAT * bufferTop,
            int lowOffset) :
        BoundaryStencil<FlowField>(parameters),
        _bufferLeft(bufferLeft), _bufferRight(bufferRight),
        _bufferBottom(bufferBottom), _bufferTop(bufferTop),
        _lowOffset(lowOffset) {}

    /** Constructor for 3D case
     * @param parameters Parameters of the simulation
     * @param bufferLeft Pointer to an array of length (N+2)^2 for saving the values at the left boundary
     * @param bufferRight Pointer to an array of length (N+2)^2 for saving the values at the right boundary
     * @param bufferBottom Pointer to an array of length (N+2)^2 for saving the values at the bottom boundary
     * @param bufferTop Pointer to an array of length (N+2)^2 for saving the values at the top boundary
     * @param bufferFront Pointer to an array of length (N+2)^2 for saving the values at the front boundary
     * @param bufferBack Pointer to an array of length (N+2)^2 for saving the values at the back boundary
     */
    ScalarBufferFillStencil(const Parameters & parameters,
            FLOAT * bufferLeft, FLOAT * bufferRight,
            FLOAT * bufferBottom, FLOAT * bufferTop,
            FLOAT * bufferFront, FLOAT * bufferBack,
            int lowOffset) :
        BoundaryStencil<FlowField>(parameters),
        _bufferLeft(bufferLeft), _bufferRight(bufferRight),
        _bufferBottom(bufferBottom), _bufferTop(bufferTop),
        _bufferFront(bufferFront), _bufferBack(bufferBack),
        _lowOffset(lowOffset) {}

    // 2D problem

    void applyLeftWall   ( FlowField & flowField, int i, int j ) {
        _bufferLeft[j - _lowOffset] = getScalar(flowField, i+2, j);
    }

    void applyRightWall  ( FlowField & flowField, int i, int j ) {
        _bufferRight[(j - _lowOffset)] = getScalar(flowField, i-1, j);
    }

    void applyBottomWall ( FlowField & flowField, int i, int j ) {
        _bufferBottom[i - _lowOffset] = getScalar(flowField, i, j+2);
    }

    void applyTopWall    ( FlowField & flowField, int i, int j ) {
        _bufferTop[(i - _lowOffset)] = getScalar(flowField, i, j-1);
    }


    // 3D problem

    void applyLeftWall   ( FlowField & flowField, int i, int j, int k ) {
        _bufferLeft[(j - _lowOffset) * (BoundaryStencil<FlowField>::_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = getScalar(flowField, i+2, j, k);
    }

    void applyRightWall  ( FlowField & flowField, int i, int j, int k ) {
        _bufferRight[(j - _lowOffset) * (BoundaryStencil<FlowField>::_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = getScalar(flowField, i-1, j, k);

    }

    void applyBottomWall ( FlowField & flowField, int i, int j, int k ) {
        _bufferBottom[(i - _lowOffset) * (BoundaryStencil<FlowField>::_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = getScalar(flowField, i, j+2, k);
    }

    void applyTopWall    ( FlowField & flowField, int i, int j, int k ) {
        _bufferTop[(i - _lowOffset) * (BoundaryStencil<FlowField>::_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = getScalar(flowField, i, j-1, k);
    }

    void applyFrontWall  ( FlowField & flowField, int i, int j, int k ) {
        _bufferFront[(i - _lowOffset) * (BoundaryStencil<FlowField>::_parameters.parallel.localSize[1] + 2) + (j - _lowOffset)] = getScalar(flowField, i, j, k+2);
    }

    void applyBackWall   ( FlowField & flowField, int i, int j, int k ) {
        _bufferBack[(i - _lowOffset) * (BoundaryStencil<FlowField>::_parameters.parallel.localSize[1] + 2) + (j - _lowOffset)] = getScalar(flowField, i, j, k-1);
    }

protected:

    virtual FLOAT & getScalar( FlowField & flowField, int i, int j ) = 0;
    virtual FLOAT & getScalar( FlowField & flowField, int i, int j, int k ) = 0;

};

#endif
