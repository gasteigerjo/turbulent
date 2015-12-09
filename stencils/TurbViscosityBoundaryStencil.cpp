#include "TurbViscosityBoundaryStencil.h"

BFInputTurbViscosityStencil::BFInputTurbViscosityStencil (const Parameters & parameters) :
    BoundaryStencil<TurbulentFlowField> (parameters),
    // Here, the obstacle size is set to zero if it was set as negative at the configuration
    _stepSize (parameters.bfStep.yRatio > 0.0 ? parameters.bfStep.yRatio*parameters.geometry.lengthY : 0.0)
{

}

// Most of the functions are empty, and they shouldn't be called, assuming that the input is always
// located at the left.

void BFInputTurbViscosityStencil::applyLeftWall   ( TurbulentFlowField & turbFlowField, int i, int j ){
    turbFlowField.getTurbViscosity().getScalar(i,j) = turbFlowField.getTurbViscosity().getScalar(i+1,j);
}

void BFInputTurbViscosityStencil::applyRightWall  ( TurbulentFlowField & turbFlowField, int i, int j ){}
void BFInputTurbViscosityStencil::applyBottomWall ( TurbulentFlowField & turbFlowField, int i, int j ){}
void BFInputTurbViscosityStencil::applyTopWall    ( TurbulentFlowField & turbFlowField, int i, int j ){}

void BFInputTurbViscosityStencil::applyLeftWall   ( TurbulentFlowField & turbFlowField, int i, int j, int k ){
    turbFlowField.getTurbViscosity().getScalar(i,j,k) = turbFlowField.getTurbViscosity().getScalar(i+1,j,k);
}

void BFInputTurbViscosityStencil::applyRightWall  ( TurbulentFlowField & turbFlowField, int i, int j, int k ){}
void BFInputTurbViscosityStencil::applyBottomWall ( TurbulentFlowField & turbFlowField, int i, int j, int k ){}
void BFInputTurbViscosityStencil::applyTopWall    ( TurbulentFlowField & turbFlowField, int i, int j, int k ){}
void BFInputTurbViscosityStencil::applyFrontWall  ( TurbulentFlowField & turbFlowField, int i, int j, int k ){}
void BFInputTurbViscosityStencil::applyBackWall   ( TurbulentFlowField & turbFlowField, int i, int j, int k ){}



NeumannTurbViscosityBoundaryStencil::NeumannTurbViscosityBoundaryStencil(const Parameters & parameters):
    BoundaryStencil<TurbulentFlowField>(parameters) {}


void NeumannTurbViscosityBoundaryStencil::applyLeftWall   ( TurbulentFlowField & turbFlowField, int i, int j ){

}

void NeumannTurbViscosityBoundaryStencil::applyRightWall  ( TurbulentFlowField & turbFlowField, int i, int j ){
    turbFlowField.getTurbViscosity().getScalar(i,j) = turbFlowField.getTurbViscosity().getScalar(i-1,j);
}

void NeumannTurbViscosityBoundaryStencil::applyBottomWall ( TurbulentFlowField & turbFlowField, int i, int j ){

}

void NeumannTurbViscosityBoundaryStencil::applyTopWall    ( TurbulentFlowField & turbFlowField, int i, int j ){

}

void NeumannTurbViscosityBoundaryStencil::applyLeftWall   ( TurbulentFlowField & turbFlowField, int i, int j, int k ){

}

void NeumannTurbViscosityBoundaryStencil::applyRightWall  ( TurbulentFlowField & turbFlowField, int i, int j, int k ){
    turbFlowField.getTurbViscosity().getScalar(i,j,k) = turbFlowField.getTurbViscosity().getScalar(i-1,j,k);
}

void NeumannTurbViscosityBoundaryStencil::applyBottomWall ( TurbulentFlowField & turbFlowField, int i, int j, int k ){

}

void NeumannTurbViscosityBoundaryStencil::applyTopWall    ( TurbulentFlowField & turbFlowField, int i, int j, int k ){

}

void NeumannTurbViscosityBoundaryStencil::applyFrontWall  ( TurbulentFlowField & turbFlowField, int i, int j, int k ){

}

void NeumannTurbViscosityBoundaryStencil::applyBackWall   ( TurbulentFlowField & turbFlowField, int i, int j, int k ){

}



MovingWallTurbViscosityStencil::MovingWallTurbViscosityStencil ( const Parameters & parameters ) :
    BoundaryStencil<TurbulentFlowField> ( parameters ) {}

// 2D stencils

void MovingWallTurbViscosityStencil::applyLeftWall ( TurbulentFlowField & turbFlowField, int i, int j ){
    turbFlowField.getTurbViscosity().getScalar(i,j) = -turbFlowField.getTurbViscosity().getScalar(i+1,j);
}


void MovingWallTurbViscosityStencil::applyRightWall ( TurbulentFlowField & turbFlowField, int i, int j ){
    turbFlowField.getTurbViscosity().getScalar(i,j) = -turbFlowField.getTurbViscosity().getScalar(i-1,j);
}


void MovingWallTurbViscosityStencil::applyBottomWall ( TurbulentFlowField & turbFlowField, int i, int j ){
    turbFlowField.getTurbViscosity().getScalar(i,j) = -turbFlowField.getTurbViscosity().getScalar(i,j+1);
}


void MovingWallTurbViscosityStencil::applyTopWall ( TurbulentFlowField & turbFlowField, int i, int j ){
    turbFlowField.getTurbViscosity().getScalar(i,j) = -turbFlowField.getTurbViscosity().getScalar(i,j-1);
}


// 3D stencils

void MovingWallTurbViscosityStencil::applyLeftWall ( TurbulentFlowField & turbFlowField, int i, int j, int k ){
    turbFlowField.getTurbViscosity().getScalar(i,j,k) = -turbFlowField.getTurbViscosity().getScalar(i+1,j,k);
}


void MovingWallTurbViscosityStencil::applyRightWall ( TurbulentFlowField & turbFlowField, int i, int j , int k ){
    turbFlowField.getTurbViscosity().getScalar(i,j,k) = -turbFlowField.getTurbViscosity().getScalar(i-1,j,k);
}


void MovingWallTurbViscosityStencil::applyBottomWall ( TurbulentFlowField & turbFlowField, int i, int j, int k ){
    turbFlowField.getTurbViscosity().getScalar(i,j,k) = -turbFlowField.getTurbViscosity().getScalar(i,j+1,k);
}


void MovingWallTurbViscosityStencil::applyTopWall ( TurbulentFlowField & turbFlowField, int i, int j, int k ){
    turbFlowField.getTurbViscosity().getScalar(i,j,k) = -turbFlowField.getTurbViscosity().getScalar(i,j-1,k);
}


void MovingWallTurbViscosityStencil::applyFrontWall ( TurbulentFlowField & turbFlowField, int i, int j, int k ){
    turbFlowField.getTurbViscosity().getScalar(i,j,k) = -turbFlowField.getTurbViscosity().getScalar(i,j,k+1);
}


void MovingWallTurbViscosityStencil::applyBackWall ( TurbulentFlowField & turbFlowField, int i, int j, int k ){
    turbFlowField.getTurbViscosity().getScalar(i,j,k) = -turbFlowField.getTurbViscosity().getScalar(i,j,k11);
}
