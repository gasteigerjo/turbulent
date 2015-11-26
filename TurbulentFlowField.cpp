#include "TurbulentFlowField.h"

TurbulentFlowField::TurbulentFlowField(const Parameters & parameters):
    FlowField(parameters),
    _turbViscosity(parameters.geometry.dim==2?ScalarField(getCellsX() + 3, getCellsY() + 3):
                      ScalarField(getCellsX() + 3, getCellsY() + 3, getCellsZ() + 3)),
    _distNearestWall(parameters.geometry.dim==2?ScalarField(getCellsX() + 3, getCellsY() + 3):
                      ScalarField(getCellsX() + 3, getCellsY() + 3, getCellsZ() + 3))
  {}

ScalarField & TurbulentFlowField::getTurbViscosity () {
    return _turbViscosity;
}

ScalarField & TurbulentFlowField::getDistNearestWall () {
    return _distNearestWall;
}
