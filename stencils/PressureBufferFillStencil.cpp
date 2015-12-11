#include "PressureBufferFillStencil.h"

PressureBufferFillStencil::PressureBufferFillStencil(const Parameters & parameters,
        FLOAT * pressuresLeft, FLOAT * pressuresRight,
        FLOAT * pressuresBottom, FLOAT * pressuresTop,
        int lowOffset) :
    ScalarBufferFillStencil<FlowField>(parameters,pressuresLeft,pressuresRight,pressuresBottom,pressuresTop,lowOffset)
    {}
PressureBufferFillStencil::PressureBufferFillStencil(const Parameters & parameters,
        FLOAT * pressuresLeft, FLOAT * pressuresRight,
        FLOAT * pressuresBottom, FLOAT * pressuresTop,
        FLOAT * pressuresFront, FLOAT * pressuresBack,
        int lowOffset) :
    ScalarBufferFillStencil<FlowField>(parameters,pressuresLeft,pressuresRight,pressuresBottom,pressuresTop,pressuresFront,pressuresBack,lowOffset)
    {}

// 2D problem
FLOAT & PressureBufferFillStencil::getScalar( FlowField & flowField, int i, int j ){
    return flowField.getPressure().getScalar(i,j);
}

// 3D problem
FLOAT & PressureBufferFillStencil::getScalar( FlowField & flowField, int i, int j, int k ){
    return flowField.getPressure().getScalar(i,j,k);
}
