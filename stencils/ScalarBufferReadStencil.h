#ifndef _SCALAR_BUFFER_READ_STENCIL_H_
#define _SCALAR_BUFFER_READ_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include <string>

/**
 *
 * Reads the values at the boundary from an array.
 */
template<class FlowField>
class ScalarBufferReadStencil : public BoundaryStencil<FlowField> {

protected:
    FLOAT * _bufferLeft;         // Array for reading the values at the left boundary
    FLOAT * _bufferRight;        // Array for reading the values at the right boundary
    FLOAT * _bufferBottom;       // Array for reading the values at the bottom boundary
    FLOAT * _bufferTop;          // Array for reading the values at the top boundary
    FLOAT * _bufferFront;        // Array for reading the values at the front boundary
    FLOAT * _bufferBack;         // Array for reading the values at the back boundary
    int _lowOffset;

public:

    /** Constructor for 2D case
     * @param parameters Parameters of the simulation
     * @param bufferLeft Pointer to an array of length (N+2) for reading the values at the left boundary
     * @param bufferRight Pointer to an array of length (N+2) for reading the values at the right boundary
     * @param bufferBottom Pointer to an array of length (N+2) for reading the values at the bottom boundary
     * @param bufferTop Pointer to an array of length (N+2) for reading the values at the top boundary
     */
    ScalarBufferReadStencil(const Parameters & parameters,
            FLOAT * bufferLeft, FLOAT * bufferRight,
            FLOAT * bufferBottom, FLOAT * bufferTop,
            int lowOffset) :
        BoundaryStencil<FlowField>(parameters),
        _bufferLeft(bufferLeft), _bufferRight(bufferRight),
        _bufferBottom(bufferBottom), _bufferTop(bufferTop),
        _lowOffset(lowOffset) {}

    /** Constructor for 3D case
     * @param parameters Parameters of the simulation
     * @param bufferLeft Pointer to an array of length (N+2)^2 for reading the values at the left boundary
     * @param bufferRight Pointer to an array of length (N+2)^2 for reading the values at the right boundary
     * @param bufferBottom Pointer to an array of length (N+2)^2 for reading the values at the bottom boundary
     * @param bufferTop Pointer to an array of length (N+2)^2 for reading the values at the top boundary
     * @param bufferFront Pointer to an array of length (N+2)^2 for reading the values at the front boundary
     * @param bufferBack Pointer to an array of length (N+2)^2 for reading the values at the back boundary
     */
    ScalarBufferReadStencil(const Parameters & parameters,
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
        getScalar(flowField, i+1, j) = _bufferLeft[2*(j - _lowOffset)+1];
    }
    void applyRightWall  ( FlowField & flowField, int i, int j ) {
        getScalar(flowField, i, j) = _bufferRight[j - _lowOffset];
    }
    void applyBottomWall ( FlowField & flowField, int i, int j ) {
        getScalar(flowField, i, j+1) = _bufferBottom[(i - _lowOffset)];
    }
    void applyTopWall    ( FlowField & flowField, int i, int j ) {
        getScalar(flowField, i, j) = _bufferTop[i - _lowOffset];
    }


    // 3D problem

    // TODO Check if array index is right
    void applyLeftWall   ( FlowField & flowField, int i, int j, int k ) {
        getScalar(flowField, i+1, j, k) = _bufferLeft[(j - _lowOffset) * (BoundaryStencil<FlowField>::_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)];
    }
    void applyRightWall  ( FlowField & flowField, int i, int j, int k ) {
        getScalar(flowField, i, j, k) = _bufferRight[(j - _lowOffset) * (BoundaryStencil<FlowField>::_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)]; // k is the inner loop
    }
    void applyBottomWall ( FlowField & flowField, int i, int j, int k ) {
        getScalar(flowField, i, j+1, k) = _bufferBottom[(i - _lowOffset) * (BoundaryStencil<FlowField>::_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)];
    }
    void applyTopWall    ( FlowField & flowField, int i, int j, int k ) {
        getScalar(flowField, i, j, k) = _bufferTop[(i - _lowOffset) * (BoundaryStencil<FlowField>::_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)]; // k is the inner loop
    }
    void applyFrontWall  ( FlowField & flowField, int i, int j, int k ) {
        getScalar(flowField, i, j, k+1) = _bufferFront[(i - _lowOffset) * (BoundaryStencil<FlowField>::_parameters.parallel.localSize[1] + 2) + (j - _lowOffset)];
    }
    void applyBackWall   ( FlowField & flowField, int i, int j, int k ) {
        getScalar(flowField, i, j, k) = _bufferBack[(i - _lowOffset) * (BoundaryStencil<FlowField>::_parameters.parallel.localSize[1] + 2) + (j - _lowOffset)]; // j is the inner loop
    }

protected:

    virtual FLOAT & getScalar( FlowField & flowField, int i, int j ) = 0;
    virtual FLOAT & getScalar( FlowField & flowField, int i, int j, int k ) = 0;

};

#endif
