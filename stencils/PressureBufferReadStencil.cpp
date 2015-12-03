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
    flowField.getPressure().getScalar(i, j) = _pressuresLeft[j - lowOffset];
}
void PressureBufferReadStencil::applyRightWall  ( FlowField & flowField, int i, int j ) {
    flowField.getPressure().getScalar(i, j) = _pressuresRight[j - lowOffset];
}
void PressureBufferReadStencil::applyBottomWall ( FlowField & flowField, int i, int j ) {
    flowField.getPressure().getScalar(i, j) = _pressuresBottom[i - lowOffset];
}
void PressureBufferReadStencil::applyTopWall    ( FlowField & flowField, int i, int j ) {
    flowField.getPressure().getScalar(i, j) = _pressuresTop[i - lowOffset];
}


// 3D problem

// TODO Check if array index is right
void PressureBufferReadStencil::applyLeftWall   ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresLeft[(j - lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - lowOffset)]; // k is the inner loop
}
void PressureBufferReadStencil::applyRightWall  ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresRight[(j - lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - lowOffset)]; // k is the inner loop
}
void PressureBufferReadStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresBottom[(i - lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - lowOffset)]; // k is the inner loop
}
void PressureBufferReadStencil::applyTopWall    ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresTop[(i - lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - lowOffset)]; // k is the inner loop
}
void PressureBufferReadStencil::applyFrontWall  ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresFront[(i - lowOffset) * (_parameters.parallel.localSize[1] + 2) + (j - lowOffset)]; // j is the inner loop
}
void PressureBufferReadStencil::applyBackWall   ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresBack[(i - lowOffset) * (_parameters.parallel.localSize[1] + 2) + (j - lowOffset)]; // j is the inner loop
}
