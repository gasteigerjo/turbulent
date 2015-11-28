#include "TurbViscosityStencil.h"
#include "StencilFunctions.h"

TurbViscosityStencil::TurbViscosityStencil ( const Parameters & parameters ) : FieldStencil<TurbulentFlowField> ( parameters ) {}


void TurbViscosityStencil::apply ( TurbulentFlowField & turbFlowField, int i, int j ){

    const int obstacle = turbFlowField.getFlags().getValue(i, j);

    if ((obstacle & OBSTACLE_SELF) == 0){ // If this is a fluid cell
        // TODO 2D implementation
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
        FLOAT mixingLength = 0.0;
        // h is the distance to the nearest wall
        const FLOAT h = turbFlowField.getDistNearestWall().getScalar(i, j, k);

        switch (_parameters.turbulenceModel.mixingLengthModel.deltaType) {
          case IgnoreDelta:
            // do not incorporate the boundary layer thickness
            mixingLength = _parameters.turbulenceModel.mixingLengthModel.kappa * h;
            break;
          default:
            mixingLength = _parameters.turbulenceModel.mixingLengthModel.kappa * h; // TODO
            break;
        }

        // calculate the turb visc from the mixing length model
        turbFlowField.getTurbViscosity().getScalar(i, j, k) = pow(mixingLength, 2) * sqrt(2.0 * SijSij);
    }
}
