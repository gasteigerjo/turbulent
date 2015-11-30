#include "MinDtStencil.h"
#include <algorithm>
#include <math.h>


MinDtStencil::MinDtStencil (const Parameters & parameters) :
    FieldStencil<TurbulentFlowField> (parameters) {
    reset();
}

void MinDtStencil::apply (TurbulentFlowField & turbulentFlowField, int i, int j){
    FLOAT localValue = 0.5 / (1.0 / _parameters.flow.Re + turbulentFlowField.getTurbViscosity().getScalar(i, j))
      / (1.0 / pow(_parameters.meshsize->getDx(i,j), 2)
       + 1.0 / pow(_parameters.meshsize->getDy(i,j), 2) );

    _minValue = std::min(_minValue, localValue);
}

void MinDtStencil::apply (TurbulentFlowField & turbulentFlowField, int i, int j, int k){

}

void MinDtStencil::reset () {
    _minValue = MY_FLOAT_MAX;
}

const FLOAT MinDtStencil::getMinValue(){
    return _minValue;
}
