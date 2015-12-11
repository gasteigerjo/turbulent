#include "TurbViscosityBufferReadStencil.h"

TurbViscosityBufferReadStencil::TurbViscosityBufferReadStencil(const Parameters & parameters,
        FLOAT * turbViscLeft, FLOAT * turbViscRight,
        FLOAT * turbViscBottom, FLOAT * turbViscTop,
        int lowOffset) :
    ScalarBufferReadStencil<TurbulentFlowField>(parameters,turbViscLeft,turbViscRight,turbViscBottom,turbViscTop,lowOffset)
    {}
TurbViscosityBufferReadStencil::TurbViscosityBufferReadStencil(const Parameters & parameters,
        FLOAT * turbViscLeft, FLOAT * turbViscRight,
        FLOAT * turbViscBottom, FLOAT * turbViscTop,
        FLOAT * turbViscFront, FLOAT * turbViscBack,
        int lowOffset) :
    ScalarBufferReadStencil<TurbulentFlowField>(parameters,turbViscLeft,turbViscRight,turbViscBottom,turbViscTop,turbViscFront,turbViscBack,lowOffset)
    {}

// 2D problem
FLOAT & TurbViscosityBufferReadStencil::getScalar( TurbulentFlowField & turbFlowField, int i, int j ){
    return turbFlowField.getTurbViscosity().getScalar(i,j);
}

// 3D problem
FLOAT & TurbViscosityBufferReadStencil::getScalar( TurbulentFlowField & turbFlowField, int i, int j, int k ){
    return turbFlowField.getTurbViscosity().getScalar(i,j,k);
}
