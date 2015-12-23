#include "TurbViscosityStencil.h"
#include "StencilFunctions.h"

TurbViscosityStencil::TurbViscosityStencil ( const Parameters & parameters ) : FieldStencil<TurbulentFlowField> ( parameters ) {}


void TurbViscosityStencil::apply ( TurbulentFlowField & turbFlowField, int i, int j ){

    const int obstacle = turbFlowField.getFlags().getValue(i, j);

    if ((obstacle & OBSTACLE_SELF) == 0){ // If this is a fluid cell

        loadLocalVelocity2D(  turbFlowField, _localVelocity, i, j);
        loadLocalMeshsize2D(_parameters, _localMeshsize, i, j);

        // do not test here for now since no other turbulence models are available
        // if (_parameters.turbulenceModel.type == "mixingLength")

        // Sij is a component of the shear strain tensor
        // SijSij is the sum of the elementwise squares
        FLOAT SijSij = pow(dudx(_localVelocity, _localMeshsize), 2)
                     + pow(dvdy(_localVelocity, _localMeshsize), 2)
                    //  + pow(dwdz(_localVelocity, _localMeshsize), 2)
                     + 0.5 * (
                        pow(dudy_cc(_localVelocity, _localMeshsize) + dvdx_cc(_localVelocity, _localMeshsize), 2)
                      // + pow(dudz_cc(_localVelocity, _localMeshsize) + dwdx_cc(_localVelocity, _localMeshsize), 2)
                      // + pow(dvdz_cc(_localVelocity, _localMeshsize) + dwdy_cc(_localVelocity, _localMeshsize), 2)
                       );

        // the mixing length of Prandtl's model
        FLOAT mixingLength = getMixingLength(turbFlowField, i, j);

        // calculate the turb visc from the mixing length model
        turbFlowField.getTurbViscosity().getScalar(i, j) = pow(mixingLength, 2) * sqrt(2.0 * SijSij);
    }
}


void TurbViscosityStencil::apply ( TurbulentFlowField & turbFlowField, int i, int j, int k ){

    const int obstacle = turbFlowField.getFlags().getValue(i, j, k);

    if ((obstacle & OBSTACLE_SELF) == 0) { // If this is a fluid cell

        loadLocalVelocity3D(  turbFlowField, _localVelocity, i, j, k);
        loadLocalMeshsize3D(_parameters, _localMeshsize, i, j, k);

        // do not test here for now since no other turbulence models are available
        // if (_parameters.turbulenceModel.type == "mixingLength")

        // Sij is a component of the shear strain tensor
        // SijSij is the sum of the elementwise squares
        FLOAT SijSij = pow(dudx(_localVelocity, _localMeshsize), 2)
                     + pow(dvdy(_localVelocity, _localMeshsize), 2)
                     + pow(dwdz(_localVelocity, _localMeshsize), 2)
                     + 0.5 * (
                        pow(dudy_cc(_localVelocity, _localMeshsize) + dvdx_cc(_localVelocity, _localMeshsize), 2)
                      + pow(dudz_cc(_localVelocity, _localMeshsize) + dwdx_cc(_localVelocity, _localMeshsize), 2)
                      + pow(dvdz_cc(_localVelocity, _localMeshsize) + dwdy_cc(_localVelocity, _localMeshsize), 2)
                       );



        // the mixing length of Prandtl's model
        FLOAT mixingLength = getMixingLength(turbFlowField, i, j, k);

        // calculate the turb visc from the mixing length model
        turbFlowField.getTurbViscosity().getScalar(i, j, k) = pow(mixingLength, 2) * sqrt(2.0 * SijSij);
    }
}



IgnoreDeltaTurbViscosityStencil::IgnoreDeltaTurbViscosityStencil ( const Parameters & parameters ) : TurbViscosityStencil(parameters) {}

FLOAT IgnoreDeltaTurbViscosityStencil::getMixingLength ( TurbulentFlowField & turbFlowField, int i, int j ){
    const FLOAT h = turbFlowField.getDistNearestWall().getScalar(i, j);
    const FLOAT delta = _parameters.turbulenceModel.mixingLengthModel.deltaValue;
    return std::min(_parameters.turbulenceModel.mixingLengthModel.kappa * h, 0.09 * delta);
}

FLOAT IgnoreDeltaTurbViscosityStencil::getMixingLength ( TurbulentFlowField & turbFlowField, int i, int j, int k ){
    const FLOAT h = turbFlowField.getDistNearestWall().getScalar(i, j, k);
    const FLOAT delta = _parameters.turbulenceModel.mixingLengthModel.deltaValue;
    return std::min(_parameters.turbulenceModel.mixingLengthModel.kappa * h, 0.09 * delta);
}



FixedDeltaTurbViscosityStencil::FixedDeltaTurbViscosityStencil ( const Parameters & parameters ) : TurbViscosityStencil(parameters) {}

FLOAT FixedDeltaTurbViscosityStencil::getMixingLength ( TurbulentFlowField & turbFlowField, int i, int j ){
    const FLOAT h = turbFlowField.getDistNearestWall().getScalar(i, j);
    return _parameters.turbulenceModel.mixingLengthModel.kappa * h;
}

FLOAT FixedDeltaTurbViscosityStencil::getMixingLength ( TurbulentFlowField & turbFlowField, int i, int j, int k ){
    const FLOAT h = turbFlowField.getDistNearestWall().getScalar(i, j, k);
    return _parameters.turbulenceModel.mixingLengthModel.kappa * h;
}



TurbFlatPlateTurbViscosityStencil::TurbFlatPlateTurbViscosityStencil ( const Parameters & parameters ) : TurbViscosityStencil(parameters) {}

FLOAT TurbFlatPlateTurbViscosityStencil::getMixingLength ( TurbulentFlowField & turbFlowField, int i, int j ){
    const FLOAT h = turbFlowField.getDistNearestWall().getScalar(i, j);
    const FLOAT posX = _parameters.meshsize->getPosX(i,j) + 0.5 * _parameters.meshsize->getDx(i,j);
    // Downstream Reynolds number
    const FLOAT Re_x = 1.0 * posX * _parameters.flow.Re;
    const FLOAT delta = 0.382 * posX / pow(Re_x, 0.2);
    return std::min(_parameters.turbulenceModel.mixingLengthModel.kappa * h, 0.09 * delta);
}

FLOAT TurbFlatPlateTurbViscosityStencil::getMixingLength ( TurbulentFlowField & turbFlowField, int i, int j, int k ){
    const FLOAT h = turbFlowField.getDistNearestWall().getScalar(i, j, k);
    const FLOAT posX = _parameters.meshsize->getPosX(i,j,k) + 0.5 * _parameters.meshsize->getDx(i,j,k);
    // Downstream Reynolds number
    const FLOAT Re_x = 1.0 * posX * _parameters.flow.Re;
    const FLOAT delta = 0.382 * posX / pow(Re_x, 0.2);
    return std::min(_parameters.turbulenceModel.mixingLengthModel.kappa * h, 0.09 * delta);
}



BlasiusLayerTurbViscosityStencil::BlasiusLayerTurbViscosityStencil ( const Parameters & parameters ) : TurbViscosityStencil(parameters) {}

FLOAT BlasiusLayerTurbViscosityStencil::getMixingLength ( TurbulentFlowField & turbFlowField, int i, int j ){
    const FLOAT h = turbFlowField.getDistNearestWall().getScalar(i, j);
    const FLOAT posX = _parameters.meshsize->getPosX(i,j) + 0.5 * _parameters.meshsize->getDx(i,j);
    // Downstream Reynolds number
    const FLOAT Re_x = 1.0 * posX * _parameters.flow.Re;
    const FLOAT delta = 4.91 * posX / sqrt(Re_x);
    return std::min(_parameters.turbulenceModel.mixingLengthModel.kappa * h, 0.09 * delta);
}

FLOAT BlasiusLayerTurbViscosityStencil::getMixingLength ( TurbulentFlowField & turbFlowField, int i, int j, int k ){
    const FLOAT h = turbFlowField.getDistNearestWall().getScalar(i, j, k);
    const FLOAT posX = _parameters.meshsize->getPosX(i,j,k) + 0.5 * _parameters.meshsize->getDx(i,j,k);
    // Downstream Reynolds number
    const FLOAT Re_x = 1.0 * posX * _parameters.flow.Re;
    const FLOAT delta = 4.91 * posX / sqrt(Re_x);
    return std::min(_parameters.turbulenceModel.mixingLengthModel.kappa * h, 0.09 * delta);
}
