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
    // printf("Pressure Fill applyLeftWall for (i,j)=(%d, %d)... (rank %d)\n", i, j, _parameters.parallel.rank);
    // _pressuresLeft[j - _lowOffset] = flowField.getPressure().getScalar(i, j);
    _pressuresLeft[j - _lowOffset] = flowField.getPressure().getScalar(i+2, j);
}
void PressureBufferFillStencil::applyRightWall  ( FlowField & flowField, int i, int j ) {
    // printf("Pressure Fill applyRightWall for (i,j)=(%d, %d)... (rank %d)\n", i-1, j, _parameters.parallel.rank);
    // _pressuresRight[j - _lowOffset] = flowField.getPressure().getScalar(i, j);
    _pressuresRight[2*(j - _lowOffset)+1] = flowField.getPressure().getScalar(i-1, j);
    _pressuresRight[2*(j - _lowOffset)] = flowField.getPressure().getScalar(i-2, j);
}
void PressureBufferFillStencil::applyBottomWall ( FlowField & flowField, int i, int j ) {
    // printf("Pressure applyBottomWall for (i,j)=(%d, %d)... (rank %d)\n", i, j, _parameters.parallel.rank);
    // _pressuresBottom[i - _lowOffset] = flowField.getPressure().getScalar(i, j);
    _pressuresBottom[i - _lowOffset] = flowField.getPressure().getScalar(i, j+2);
}
void PressureBufferFillStencil::applyTopWall    ( FlowField & flowField, int i, int j ) {
    // printf("Pressure Fill applyTopWall for (i,j)=(%d, %d)... (rank %d)\n", i, j-1, _parameters.parallel.rank);
    // _pressuresTop[i - _lowOffset] = flowField.getPressure().getScalar(i, j);
    _pressuresTop[2*(i - _lowOffset)+1] = flowField.getPressure().getScalar(i, j-1);
    _pressuresTop[2*(i - _lowOffset)] = flowField.getPressure().getScalar(i, j-2);
}


// 3D problem

// TODO Check if array index is right
void PressureBufferFillStencil::applyLeftWall   ( FlowField & flowField, int i, int j, int k ) {
    // _pressuresLeft[(j - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = flowField.getPressure().getScalar(i, j, k); // k is the inner loop
    _pressuresLeft[(j - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = flowField.getPressure().getScalar(i+2, j, k);
}
void PressureBufferFillStencil::applyRightWall  ( FlowField & flowField, int i, int j, int k ) {
    // _pressuresRight[(j - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = flowField.getPressure().getScalar(i, j, k); // k is the inner loop
    _pressuresRight[2*(j - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset) +1] = flowField.getPressure().getScalar(i-1, j, k);
    _pressuresRight[2*(j - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = flowField.getPressure().getScalar(i-2, j, k);

}
void PressureBufferFillStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ) {
    // _pressuresBottom[(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = flowField.getPressure().getScalar(i, j, k); // k is the inner loop
    _pressuresBottom[(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = flowField.getPressure().getScalar(i, j+2, k);
}
void PressureBufferFillStencil::applyTopWall    ( FlowField & flowField, int i, int j, int k ) {
    // _pressuresTop[(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = flowField.getPressure().getScalar(i, j, k); // k is the inner loop
    _pressuresTop[2*(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset) +1] = flowField.getPressure().getScalar(i, j-1, k);
    _pressuresTop[2*(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (k - _lowOffset)] = flowField.getPressure().getScalar(i, j-2, k);
}
void PressureBufferFillStencil::applyFrontWall  ( FlowField & flowField, int i, int j, int k ) {
    // _pressuresFront[(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (j - _lowOffset)] = flowField.getPressure().getScalar(i, j, k); // j is the inner loop
    _pressuresFront[(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (j - _lowOffset)] = flowField.getPressure().getScalar(i, j, k+2);
}
void PressureBufferFillStencil::applyBackWall   ( FlowField & flowField, int i, int j, int k ) {
    // _pressuresBack[(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (j - _lowOffset)] = flowField.getPressure().getScalar(i, j, k); // j is the inner loop
    _pressuresBack[2*(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (j - _lowOffset) +1] = flowField.getPressure().getScalar(i, j, k-1);
    _pressuresBack[2*(i - _lowOffset) * (_parameters.parallel.localSize[2] + 2) + (j - _lowOffset)] = flowField.getPressure().getScalar(i, j, k-2);
}
