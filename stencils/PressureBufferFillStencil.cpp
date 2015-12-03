#include "PressureBufferFillStencil.h"

PressureBufferFillStencil::PressureBufferFillStencil(const Parameters & parameters,
        FLOAT * pressuresLeft, FLOAT * pressuresRight,
        FLOAT * pressuresBottom, FLOAT * pressuresTop,
        int lowOffset) :
    BoundaryStencil<FlowField>(parameters),
    _pressuresLeft(pressuresLeft), _pressuresRight(pressuresRight),
    _pressuresBottom(pressuresBottom), _pressuresTop(pressuresTop),
    _lowOffset(lowOffset) {}

PressureBufferFillStencil::PressureBufferFillStencil(const Parameters & parameters,
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

void PressureBufferFillStencil::applyLeftWall   ( FlowField & flowField, int i, int j ) {
    printf("Pressure applyLeftWall for (i,j)=(%d, %d)... (rank %d)\n", i, j, _parameters.parallel.rank);
    _pressuresLeft[j - _lowOffset] = flowField.getPressure().getScalar(i, j);
}
void PressureBufferFillStencil::applyRightWall  ( FlowField & flowField, int i, int j ) {
    printf("Pressure applyRightWall for (i,j)=(%d, %d)... (rank %d)\n", i, j, _parameters.parallel.rank);
    _pressuresRight[j - _lowOffset] = flowField.getPressure().getScalar(i, j);
}
void PressureBufferFillStencil::applyBottomWall ( FlowField & flowField, int i, int j ) {
    printf("Pressure applyBottomWall for (i,j)=(%d, %d)... (rank %d)\n", i, j, _parameters.parallel.rank);
    _pressuresBottom[i - _lowOffset] = flowField.getPressure().getScalar(i, j);
}
void PressureBufferFillStencil::applyTopWall    ( FlowField & flowField, int i, int j ) {
    printf("Pressure applyTopWall for (i,j)=(%d, %d)... (rank %d)\n", i, j, _parameters.parallel.rank);
    _pressuresTop[i - _lowOffset] = flowField.getPressure().getScalar(i, j);
}


// 3D problem

// TODO Check if array index is right
void PressureBufferFillStencil::applyLeftWall   ( FlowField & flowField, int i, int j, int k ) {
    _pressuresLeft[(j - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = flowField.getPressure().getScalar(i, j, k); // k is the inner loop
}
void PressureBufferFillStencil::applyRightWall  ( FlowField & flowField, int i, int j, int k ) {
    _pressuresRight[(j - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = flowField.getPressure().getScalar(i, j, k); // k is the inner loop
}
void PressureBufferFillStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ) {
    _pressuresBottom[(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = flowField.getPressure().getScalar(i, j, k); // k is the inner loop
}
void PressureBufferFillStencil::applyTopWall    ( FlowField & flowField, int i, int j, int k ) {
    _pressuresTop[(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = flowField.getPressure().getScalar(i, j, k); // k is the inner loop
}
void PressureBufferFillStencil::applyFrontWall  ( FlowField & flowField, int i, int j, int k ) {
    _pressuresFront[(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (j - _lowOffset)] = flowField.getPressure().getScalar(i, j, k); // j is the inner loop
}
void PressureBufferFillStencil::applyBackWall   ( FlowField & flowField, int i, int j, int k ) {
    _pressuresBack[(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (j - _lowOffset)] = flowField.getPressure().getScalar(i, j, k); // j is the inner loop
}
