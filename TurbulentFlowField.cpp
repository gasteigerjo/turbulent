#include "TurbulentFlowField.h"

ScalarField & TurbulentFlowField::getTurbViscosity () {
    return _turbViscosity;
}

ScalarField & TurbulentFlowField::getDistNearestWall () {
    return _distNearestWall;
}
