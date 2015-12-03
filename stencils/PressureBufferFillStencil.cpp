#include "PressureBufferFillStencil.h"

PressureBufferFillStencil::PressureBufferFillStencil(const Parameters & parameters,
        FLOAT * pressuresLeft, FLOAT * pressuresRight,
        FLOAT * pressuresBottom, FLOAT * pressuresTop) :
    BoundaryStencil<FlowField>(parameters),
    _pressuresLeft(pressuresLeft), _pressuresRight(pressuresRight),
    _pressuresBottom(pressuresBottom), _pressuresTop(pressuresTop) {}

PressureBufferFillStencil::PressureBufferFillStencil(const Parameters & parameters,
        FLOAT * pressuresLeft, FLOAT * pressuresRight,
        FLOAT * pressuresBottom, FLOAT * pressuresTop,
        FLOAT * pressuresFront, FLOAT * pressuresBack) :
    BoundaryStencil<FlowField>(parameters),
    _pressuresLeft(pressuresLeft), _pressuresRight(pressuresRight),
    _pressuresBottom(pressuresBottom), _pressuresTop(pressuresTop),
    _pressuresFront(pressuresFront), _pressuresBack(pressuresBack) {}


// 2D problem

void PressureBufferFillStencil::applyLeftWall   ( FlowField & flowField, int i, int j ) {
    printf("Pressure applyLeftWall for (i,j)=(%d, %d)... (rank %d)\n", i, j, _parameters.parallel.rank);
    _pressuresLeft[j] = flowField.getPressure().getScalar(i, j);
}
void PressureBufferFillStencil::applyRightWall  ( FlowField & flowField, int i, int j ) {
    printf("Pressure applyRightWall for (i,j)=(%d, %d)... (rank %d)\n", i, j, _parameters.parallel.rank);
    _pressuresRight[j] = flowField.getPressure().getScalar(i, j);
}
void PressureBufferFillStencil::applyBottomWall ( FlowField & flowField, int i, int j ) {
    printf("Pressure applyBottomWall for (i,j)=(%d, %d)... (rank %d)\n", i, j, _parameters.parallel.rank);
    _pressuresBottom[i] = flowField.getPressure().getScalar(i, j);
}
void PressureBufferFillStencil::applyTopWall    ( FlowField & flowField, int i, int j ) {
    printf("Pressure applyTopWall for (i,j)=(%d, %d)... (rank %d)\n", i, j, _parameters.parallel.rank);
    _pressuresTop[i] = flowField.getPressure().getScalar(i, j);
}


// 3D problem

// TODO Check if array index is right
void PressureBufferFillStencil::applyLeftWall   ( FlowField & flowField, int i, int j, int k ) {
    _pressuresLeft[j * (_parameters.parallel.localSize[2] + 2) + k] = flowField.getPressure().getScalar(i, j, k); // k is the inner loop
}
void PressureBufferFillStencil::applyRightWall  ( FlowField & flowField, int i, int j, int k ) {
    _pressuresRight[j * (_parameters.parallel.localSize[2] + 2) + k] = flowField.getPressure().getScalar(i, j, k); // k is the inner loop
}
void PressureBufferFillStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ) {
    _pressuresBottom[i * (_parameters.parallel.localSize[2] + 2) + k] = flowField.getPressure().getScalar(i, j, k); // k is the inner loop
}
void PressureBufferFillStencil::applyTopWall    ( FlowField & flowField, int i, int j, int k ) {
    _pressuresTop[i * (_parameters.parallel.localSize[2] + 2) + k] = flowField.getPressure().getScalar(i, j, k); // k is the inner loop
}
void PressureBufferFillStencil::applyFrontWall  ( FlowField & flowField, int i, int j, int k ) {
    _pressuresFront[i * (_parameters.parallel.localSize[1] + 2) + j] = flowField.getPressure().getScalar(i, j, k); // j is the inner loop
}
void PressureBufferFillStencil::applyBackWall   ( FlowField & flowField, int i, int j, int k ) {
    _pressuresBack[i * (_parameters.parallel.localSize[1] + 2) + j] = flowField.getPressure().getScalar(i, j, k); // j is the inner loop
}
