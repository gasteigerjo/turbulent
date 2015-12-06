#include "PressureBufferReadStencil.h"

PressureBufferReadStencil::PressureBufferReadStencil(const Parameters & parameters,
        FLOAT * pressuresLeft, FLOAT * pressuresRight,
        FLOAT * pressuresBottom, FLOAT * pressuresTop,
        int lowOffset) :
    BoundaryStencil<FlowField>(parameters),
    _pressuresLeft(pressuresLeft), _pressuresRight(pressuresRight),
    _pressuresBottom(pressuresBottom), _pressuresTop(pressuresTop),
    _lowOffset(lowOffset) {}

PressureBufferReadStencil::PressureBufferReadStencil(const Parameters & parameters,
        FLOAT * pressuresLeft, FLOAT * pressuresRight,
        FLOAT * pressuresBottom, FLOAT * pressuresTop,
        FLOAT * pressuresFront, FLOAT * pressuresBack,
        int lowOffset) :
    BoundaryStencil<FlowField>(parameters),
    _pressuresLeft(pressuresLeft), _pressuresRight(pressuresRight),
    _pressuresBottom(pressuresBottom), _pressuresTop(pressuresTop),
    _pressuresFront(pressuresFront), _pressuresBack(pressuresBack),
    _lowOffset(lowOffset) {}


// 2D problem

void PressureBufferReadStencil::applyLeftWall   ( FlowField & flowField, int i, int j ) {
    // printf("Pressure Read applyLeftWall for (i,j)=(%d, %d)... (rank %d)\n", i, j, _parameters.parallel.rank);
    // flowField.getPressure().getScalar(i-1, j) = _pressuresLeft[j - _lowOffset];
    flowField.getPressure().getScalar(i+1, j) = _pressuresLeft[2*(j - _lowOffset)+1];
    flowField.getPressure().getScalar(i, j) = _pressuresLeft[2*(j - _lowOffset)];
}
void PressureBufferReadStencil::applyRightWall  ( FlowField & flowField, int i, int j ) {
    // flowField.getPressure().getScalar(i+1, j) = _pressuresRight[j - _lowOffset];
    flowField.getPressure().getScalar(i, j) = _pressuresRight[j - _lowOffset];
}
void PressureBufferReadStencil::applyBottomWall ( FlowField & flowField, int i, int j ) {
    // printf("Pressure Read applyBottomWall for (i,j)=(%d, %d)... (rank %d)\n", i, j, _parameters.parallel.rank);
    // flowField.getPressure().getScalar(i, j-1) = _pressuresBottom[i - _lowOffset];
    flowField.getPressure().getScalar(i, j+1) = _pressuresBottom[2*(i - _lowOffset)+1];
    flowField.getPressure().getScalar(i, j) = _pressuresBottom[2*(i - _lowOffset)];
}
void PressureBufferReadStencil::applyTopWall    ( FlowField & flowField, int i, int j ) {
    // flowField.getPressure().getScalar(i, j+1) = _pressuresTop[i - _lowOffset];
    flowField.getPressure().getScalar(i, j) = _pressuresTop[i - _lowOffset];
}


// 3D problem

// TODO Check if array index is right
void PressureBufferReadStencil::applyLeftWall   ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresLeft[(j - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)]; // k is the inner loop
}
void PressureBufferReadStencil::applyRightWall  ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresRight[(j - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)]; // k is the inner loop
}
void PressureBufferReadStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresBottom[(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)]; // k is the inner loop
}
void PressureBufferReadStencil::applyTopWall    ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresTop[(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)]; // k is the inner loop
}
void PressureBufferReadStencil::applyFrontWall  ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresFront[(i - _lowOffset) * (_parameters.parallel.localSize[1] + 2) + (j - _lowOffset)]; // j is the inner loop
}
void PressureBufferReadStencil::applyBackWall   ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresBack[(i - _lowOffset) * (_parameters.parallel.localSize[1] + 2) + (j - _lowOffset)]; // j is the inner loop
}
