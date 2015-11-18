#include "PressureBufferReadStencil.h"

PressureBufferReadStencil::PressureBufferReadStencil(const Parameters & parameters,
        FLOAT * pressuresLeft, FLOAT * pressuresRight,
        FLOAT * pressuresBottom, FLOAT * pressuresTop) :
    BoundaryStencil<FlowField>(parameters),
    _pressuresLeft(pressuresLeft), _pressuresRight(pressuresRight),
    _pressuresBottom(pressuresBottom), _pressuresTop(pressuresTop) {}

PressureBufferReadStencil::PressureBufferReadStencil(const Parameters & parameters,
        FLOAT * pressuresLeft, FLOAT * pressuresRight,
        FLOAT * pressuresBottom, FLOAT * pressuresTop,
        FLOAT * pressuresFront, FLOAT * pressuresBack) :
    BoundaryStencil<FlowField>(parameters),
    _pressuresLeft(pressuresLeft), _pressuresRight(pressuresRight),
    _pressuresBottom(pressuresBottom), _pressuresTop(pressuresTop),
    _pressuresFront(pressuresFront), _pressuresBack(pressuresBack) {}


// 2D problem

void PressureBufferReadStencil::applyLeftWall   ( FlowField & flowField, int i, int j ) {
    flowField.getPressure().getScalar(i, j) = _pressuresLeft[j];
}
void PressureBufferReadStencil::applyRightWall  ( FlowField & flowField, int i, int j ) {
    flowField.getPressure().getScalar(i, j) = _pressuresRight[j];
}
void PressureBufferReadStencil::applyBottomWall ( FlowField & flowField, int i, int j ) {
    flowField.getPressure().getScalar(i, j) = _pressuresBottom[i];
}
void PressureBufferReadStencil::applyTopWall    ( FlowField & flowField, int i, int j ) {
    flowField.getPressure().getScalar(i, j) = _pressuresTop[i];
}


// 3D problem

// TODO Check if array index is right
void PressureBufferReadStencil::applyLeftWall   ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresLeft[j * (_parameters.parallel.localSize[2] + 2) + k]; // k is the inner loop
}
void PressureBufferReadStencil::applyRightWall  ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresRight[j * (_parameters.parallel.localSize[2] + 2) + k]; // k is the inner loop
}
void PressureBufferReadStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresBottom[i * (_parameters.parallel.localSize[2] + 2) + k]; // k is the inner loop
}
void PressureBufferReadStencil::applyTopWall    ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresTop[i * (_parameters.parallel.localSize[2] + 2) + k]; // k is the inner loop
}
void PressureBufferReadStencil::applyFrontWall  ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresFront[i * (_parameters.parallel.localSize[1] + 2) + j]; // j is the inner loop
}
void PressureBufferReadStencil::applyBackWall   ( FlowField & flowField, int i, int j, int k ) {
    flowField.getPressure().getScalar(i, j, k) = _pressuresBack[i * (_parameters.parallel.localSize[1] + 2) + j]; // j is the inner loop
}
