#include "PostStencil.h"
#include "../Iterators.h"
#include <string>
#include <ostream>
#include <sstream>


BasicPostStencil::BasicPostStencil ( const Parameters & parameters )
  : PostStencil<FlowField>(parameters) {}

void BasicPostStencil::apply ( FlowField & flowField, int i, int j ) {
    FLOAT pressure = 0.0;
    FLOAT* velocity = new FLOAT(2);
    const int obstacle = flowField.getFlags().getValue(i, j);

    // check whether current cell is a fluid cell or not
    if ((obstacle & OBSTACLE_SELF) == 0) {
      flowField.getPressureAndVelocity(pressure, velocity,  i, j);
    } else {
      velocity[0] = 0.0;
      velocity[1] = 0.0;
    }

    _ssPressure << pressure << std::endl;
    _ssVelocity << velocity[0] << " " << velocity[1] << " 0" << std::endl;
    // _ssFlags << obstacle << std::endl;
}

void BasicPostStencil::apply ( FlowField & flowField, int i, int j, int k ) {
    FLOAT pressure = 0.0;
    FLOAT* velocity = new FLOAT(3);
    const int obstacle = flowField.getFlags().getValue(i, j, k);

    // check whether current cell is a fluid cell or not
    if ((obstacle & OBSTACLE_SELF) == 0) {
      flowField.getPressureAndVelocity(pressure, velocity,  i, j, k);
    } else {
      velocity[0] = 0.0;
      velocity[1] = 0.0;
      velocity[2] = 0.0;
    }

    _ssPressure << pressure << std::endl;
    _ssVelocity << velocity[0] << " " << velocity[1] << " " << velocity[2] << std::endl;
    // _ssFlags << obstacle << std::endl;
}

void BasicPostStencil::writeVtk ( std::ostream & stream ) {
    // write pressure data
    stream << "\nSCALARS pressure float 1" << std::endl;
    stream << "LOOKUP_TABLE default" << std::endl;
    stream << _ssPressure.str();

    // write velocity data
    stream << "\nVECTORS velocity float" << std::endl;
    stream << _ssVelocity.str();

    // // write flags data
    // stream << "\nSCALARS flags integer 1" << std::endl;
    // stream << "LOOKUP_TABLE default" << std::endl;
    // stream << _ssFlags.str();

    // clear the stringstreams
    _ssPressure.str("");
    _ssVelocity.str("");
    // _ssFlags.str("");
}



WallPostStencil::WallPostStencil ( const Parameters & parameters )
  : PostStencil<FlowField>(parameters),
    _sizeX(parameters.geometry.lengthX),
    _sizeY(parameters.geometry.lengthY),
    _sizeZ(parameters.geometry.lengthZ),
    _stepX(_parameters.bfStep.xRatio * parameters.geometry.lengthX),
    _stepY(_parameters.bfStep.yRatio * parameters.geometry.lengthY) { // not working if the domain is split in y and/or z direction

    // find the j-index of the first fluid cells on top of the step
    _jOnStep = -1;
    if (_parameters.meshsize->getPosY(0,1,0) < _stepY){
      // only get the index of the cells on top of the step if the step is possibly in this domain
      for (int j = 2; j < 2 + _parameters.parallel.localSize[1]; j++) {
        if (_stepY < _parameters.meshsize->getPosY(0,j,0) + 0.5 * _parameters.meshsize->getDy(0,j,0)){
          _jOnStep = j;
          break;
        }
      }
    }

    // find the i-index of the first fluid cells after the step
    _iBehindStep = -1; // indication that the step is not in this domain
    if (_parameters.meshsize->getPosX(1,0,0) < _stepX){
      // only get the index of the cells behind the step if the step is possibly in this domain
      for (int i = 2; i < 2 + _parameters.parallel.localSize[1]; i++) {
        if (_stepX < _parameters.meshsize->getPosX(i,0,0) + 0.5 * _parameters.meshsize->getDx(i,0,0)){
          _iBehindStep = i;
          break;
        }
      }
    }
}

WallPostStencil::~WallPostStencil(){
}

void WallPostStencil::preapply ( FlowField & flowField ) {

}

void WallPostStencil::apply ( FlowField & flowField, int i, int j ) {
  // not yet implemented for 2D
  // TODO: is it necessary for 2D? will we need it?
}

void WallPostStencil::apply ( FlowField & flowField, int i, int j, int k ) {
    Meshsize *ms = _parameters.meshsize;
    FLOAT tauw[3];
    tauw[0] = 0.0;
    tauw[1] = 0.0;
    tauw[2] = 0.0;

    const int obstacle = flowField.getFlags().getValue(i,j,k);
    const int globalJ = _parameters.parallel.firstCorner[1] + j - 2;
    const int globalK = _parameters.parallel.firstCorner[2] + k - 2;

    // check whether current cell is a fluid cell or not
    if ((obstacle & OBSTACLE_SELF) == 0) {
      FLOAT velocity[3];
      flowField.getVelocityCenter(velocity, i, j, k);

      const FLOAT posX = ms->getPosX(i,j,k) + 0.5*ms->getDx(i,j,k);
      const FLOAT posY = ms->getPosY(i,j,k) + 0.5*ms->getDy(i,j,k);
      // const FLOAT posZ = ms->getPosZ(i,j,k) + 0.5*ms->getDz(i,j,k);

      if ((globalJ == 0 && posX > _stepX) || (j == _jOnStep && posX < _stepX)
        || (globalJ == _parameters.geometry.sizeY - 1)) {
        // bottom/top cell
        tauw[0] = velocity[0] / (0.5 * ms->getDy(i,j,k));
        tauw[2] = velocity[2] / (0.5 * ms->getDy(i,j,k));
      } else if ((globalK == 0)
        || (globalK == _parameters.geometry.sizeZ - 1)) {
        // front/back cell
        tauw[0] = velocity[0] / (0.5 * ms->getDz(i,j,k));
        tauw[1] = velocity[1] / (0.5 * ms->getDz(i,j,k));
      } else if (i == _iBehindStep && posY < _stepY) {
        // cell behind step
        tauw[1] = velocity[1] / (0.5 * ms->getDx(i,j,k));
        tauw[2] = velocity[2] / (0.5 * ms->getDx(i,j,k));
      }
    } else {
    }

    _ssTauw << tauw[0] << " " << tauw[1] << " " << tauw[2] << std::endl;
}

void WallPostStencil::writeVtk ( std::ostream & stream ) {
    // write wall shear stress
    stream << "\nVECTORS wallShearStress float" << std::endl;
    stream << _ssTauw.str();

    // clear the stringstreams
    _ssTauw.str("");
}
