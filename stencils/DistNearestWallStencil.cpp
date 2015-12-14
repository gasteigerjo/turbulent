#include "DistNearestWallStencil.h"

DistNearestWallStencil::DistNearestWallStencil ( const Parameters & parameters ) :
    FieldStencil<TurbulentFlowField> ( parameters ),
    sizeY(parameters.geometry.lengthY),
    sizeZ(parameters.geometry.lengthZ),
    stepX(parameters.bfStep.xRatio*parameters.geometry.lengthX),
    stepY(parameters.bfStep.yRatio*parameters.geometry.lengthY)
{}


void DistNearestWallStencil::apply ( TurbulentFlowField & turbFlowField, int i, int j ){

    const int obstacle = turbFlowField.getFlags().getValue(i, j);

    if ((obstacle & OBSTACLE_SELF) == 0){ // If this is a fluid cell
        // position of the cell center
        const FLOAT posX = _parameters.meshsize->getPosX(i,j) + 0.5 * _parameters.meshsize->getDx(i,j);
        const FLOAT posY = _parameters.meshsize->getPosY(i,j) + 0.5 * _parameters.meshsize->getDy(i,j);

        FLOAT distance = 0.0;

        // calculate the NORMAL distance to the nearest wall
        if (_parameters.simulation.scenario == "channel"){
          // check where the cell is located relative to the step
          if (posX <= stepX){
            // cell is located above the step
            distance = std::min(sizeY - posY, posY - stepY);
          }else if (posX < stepX + stepY && posY <= stepY){
            // cell is located behind the step in its proximity
            // decide whether nearest wall is top/bottom wall or vertical wall of the step
            distance = std::min(sizeY - posY, posY);
            distance = std::min(distance, posX - stepX);
          }else{
            // cell is located well downstream of the step and
            // the distance is thus not influenced by the step
            distance = std::min(sizeY - posY, posY);
          }
        }else{
          handleError(1, "Only channel scenario supported for turbulence simulation");
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

        // calculate the NORMAL distance to the nearest wall
        if (_parameters.simulation.scenario == "channel"){
          // check where the cell is located relative to the step
          if (posX <= stepX){
            // cell is located above the step
            distance = std::min(sizeY - posY, posY - stepY);
          }else if (posX < stepX + stepY && posY <= stepY){
            // cell is located behind the step in its proximity
            // decide whether nearest wall is top/bottom wall or vertical wall of the step
            distance = std::min(sizeY - posY, posY);
            distance = std::min(distance, posX - stepX);
          }else{
            // cell is located well downstream of the step and
            // the distance is thus not influenced by the step
            distance = std::min(sizeY - posY, posY);
          }

          distance = std::min(distance, sizeZ - posZ);
          distance = std::min(distance, posZ);
        }else{
          handleError(1, "Only channel scenario supported for turbulence simulation");
        }

        turbFlowField.getDistNearestWall().getScalar(i, j, k) = distance;
    }
}
