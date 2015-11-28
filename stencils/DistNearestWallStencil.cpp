#include "DistNearestWallStencil.h"

DistNearestWallStencil::DistNearestWallStencil ( const Parameters & parameters ) :
    FieldStencil<TurbulentFlowField> ( parameters ),
    sizeY(parameters.geometry.lengthY),
    sizeZ(parameters.geometry.lengthZ),
    stepX(parameters.bfStep.xRatio*parameters.geometry.lengthX),
    stepY(parameters.bfStep.xRatio*parameters.geometry.lengthY)
{}


void DistNearestWallStencil::apply ( TurbulentFlowField & turbFlowField, int i, int j ){

    const int obstacle = turbFlowField.getFlags().getValue(i, j);

    if ((obstacle & OBSTACLE_SELF) == 0){ // If this is a fluid cell
        // position of the cell center
        const FLOAT posX = _parameters.meshsize->getPosX(i,j) + 0.5 * _parameters.meshsize->getDx(i,j);
        const FLOAT posY = _parameters.meshsize->getPosY(i,j) + 0.5 * _parameters.meshsize->getDy(i,j);

        FLOAT distance = 0.0;

        if (_parameters.simulation.scenario == "channel"){
          // check whether cell is located above the step
          if (posX <= stepX){
            distance = std::min(sizeY - posY, posY - stepY);
          }else{
            distance = std::min(sizeY - posY, posY);
          }
        }else{
          // TODO what to do in other scenarios?
        }

        turbFlowField.getDistNearestWall().getScalar(i, j) = distance;
    }
}


void DistNearestWallStencil::apply ( TurbulentFlowField & turbFlowField, int i, int j, int k ){

    const int obstacle = turbFlowField.getFlags().getValue(i, j, k);

    if ((obstacle & OBSTACLE_SELF) == 0) { // If this is a fluid cell
        // position of the cell center
        const FLOAT posX = _parameters.meshsize->getPosX(i,j,k) + 0.5 * _parameters.meshsize->getDx(i,j,k);
        const FLOAT posY = _parameters.meshsize->getPosY(i,j,k) + 0.5 * _parameters.meshsize->getDy(i,j,k);
        const FLOAT posZ = _parameters.meshsize->getPosZ(i,j,k) + 0.5 * _parameters.meshsize->getDz(i,j,k);

        FLOAT distance = 0.0;

        if (_parameters.simulation.scenario == "channel"){
          // check whether cell is located above the step
          if (posX <= stepX){
            distance = std::min(sizeY - posY, posY - stepY);
          }else{
            distance = std::min(sizeY - posY, posY);
          }

          distance = std::min(distance, sizeZ - posZ);
          distance = std::min(distance, posZ);
        }else{
          // TODO what to do in other scenarios?
        }

        turbFlowField.getDistNearestWall().getScalar(i, j, k) = distance;
    }
}
