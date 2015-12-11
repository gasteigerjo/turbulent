#include "PressureBufferReadStencil.h"

PressureBufferReadStencil::PressureBufferReadStencil(const Parameters & parameters,
        FLOAT * pressuresLeft, FLOAT * pressuresRight,
        FLOAT * pressuresBottom, FLOAT * pressuresTop,
        int lowOffset) :
    ScalarBufferReadStencil<FlowField>(parameters,pressuresLeft,pressuresRight,pressuresBottom,pressuresTop,lowOffset)
    {}
PressureBufferReadStencil::PressureBufferReadStencil(const Parameters & parameters,
        FLOAT * pressuresLeft, FLOAT * pressuresRight,
        FLOAT * pressuresBottom, FLOAT * pressuresTop,
        FLOAT * pressuresFront, FLOAT * pressuresBack,
        int lowOffset) :
    ScalarBufferReadStencil<FlowField>(parameters,pressuresLeft,pressuresRight,pressuresBottom,pressuresTop,pressuresFront,pressuresBack,lowOffset)
    {}

// 2D problem
FLOAT & PressureBufferReadStencil::getScalar( FlowField & flowField, int i, int j ){
    return flowField.getPressure().getScalar(i,j);
}

// 3D problem
FLOAT & PressureBufferReadStencil::getScalar( FlowField & flowField, int i, int j, int k ){
    return flowField.getPressure().getScalar(i,j,k);
}
