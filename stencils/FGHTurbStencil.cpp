#include "FGHTurbStencil.h"
#include "StencilFunctions.h"
#include "Definitions.h"

FGHTurbStencil::FGHTurbStencil ( const Parameters & parameters ) : FieldStencil<TurbulentFlowField> ( parameters ) {}


void FGHTurbStencil::apply ( TurbulentFlowField & turbulentFlowField,  int i, int j ){

    // Load local velocities into the center layer of the local array

    loadLocalVelocity2D(  turbulentFlowField, _localVelocity, i, j);
    loadLocalMeshsize2D(_parameters, _localMeshsize, i, j);
    loadLocalTurbViscosity2D(turbulentFlowField, _localTurbViscosity,i,j);

    FLOAT* const values = turbulentFlowField.getFGH().getVector(i,j);

    // Now the localVelocity array should contain lexicographically ordered elements around the
    // given index 2D

    values [0] = computeTurbF2D(_localVelocity, _localTurbViscosity, _localMeshsize, _parameters, _parameters.timestep.dt);
    values [1] = computeTurbG2D(_localVelocity, _localTurbViscosity, _localMeshsize, _parameters, _parameters.timestep.dt);

}


void FGHTurbStencil::apply ( TurbulentFlowField & turbulentFlowField, int i, int j, int k ){
    // The same as in 2D, with slight modifications

    const int obstacle = turbulentFlowField.getFlags().getValue(i, j, k);

    FLOAT * const values = turbulentFlowField.getFGH().getVector(i,j,k);

    if ((obstacle & OBSTACLE_SELF) == 0){   // If the cell is fluid

        loadLocalVelocity3D(  turbulentFlowField, _localVelocity, i, j, k);
        loadLocalMeshsize3D(_parameters, _localMeshsize, i, j, k);
  		loadLocalTurbViscosity3D( turbulentFlowField, _localTurbViscosity,i,j,k);

        if ((obstacle & OBSTACLE_RIGHT) == 0) { // If the right cell is fluid
            values [0] = computeTurbF3D(_localVelocity, _localTurbViscosity, _localMeshsize, _parameters, _parameters.timestep.dt);
        }
        if ((obstacle & OBSTACLE_TOP) == 0) {
            values [1] = computeTurbG3D(_localVelocity, _localTurbViscosity, _localMeshsize, _parameters, _parameters.timestep.dt);
        }
        if ((obstacle & OBSTACLE_BACK) == 0) {
            values [2] = computeTurbH3D(_localVelocity, _localTurbViscosity, _localMeshsize, _parameters, _parameters.timestep.dt);
        }
    }
}
