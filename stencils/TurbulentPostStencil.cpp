#include "TurbulentPostStencil.h"
#include "../Iterators.h"
#include <string>
#include <ostream>
#include <sstream>


void TurbulentPostStencil::apply ( TurbulentFlowField & turbFlowField, int i, int j ) {
    FLOAT turbVisc = 0.0;
    // FLOAT distWall = 0.0;
    const int obstacle = turbFlowField.getFlags().getValue(i, j);

    // check whether current cell is a fluid cell or not
    if ((obstacle & OBSTACLE_SELF) == 0) {
        turbVisc = turbFlowField.getTurbViscosity().getScalar(i, j);
        // distWall = turbFlowField.getDistNearestWall().getScalar(i, j);
    }
    _ssTurbVisc << turbVisc << std::endl;
    // _ssDistWall << distWall << std::endl;
}

void TurbulentPostStencil::apply ( TurbulentFlowField & turbFlowField, int i, int j, int k ) {
    FLOAT turbVisc = 0.0;
    // FLOAT distWall = 0.0;
    const int obstacle = turbFlowField.getFlags().getValue(i, j, k);

    // check whether current cell is a fluid cell or not
    if ((obstacle & OBSTACLE_SELF) == 0) {
        turbVisc = turbFlowField.getTurbViscosity().getScalar(i, j, k);
        // distWall = turbFlowField.getDistNearestWall().getScalar(i, j, k);
    }
    _ssTurbVisc << turbVisc << std::endl;
    // _ssDistWall << distWall << std::endl;
}

void TurbulentPostStencil::writeVtk ( std::ostream & stream ) {
    // write turbulent viscosity data
    stream << "\nSCALARS turbViscosity float 1" << std::endl;
    stream << "LOOKUP_TABLE default" << std::endl;
    stream << _ssTurbVisc.str();

    // // write wall distance data
    // stream << "\nSCALARS distWall float 1" << std::endl;
    // stream << "LOOKUP_TABLE default" << std::endl;
    // stream << _ssDistWall.str();

    // clear the stringstreams
    _ssTurbVisc.str("");
    // _ssDistWall.str("");
}


TurbulentWallPostStencil::TurbulentWallPostStencil ( const Parameters & parameters )
  : PostStencil<TurbulentFlowField>(parameters),
    _sizeX(parameters.geometry.lengthX),
    _sizeY(parameters.geometry.lengthY),
    _sizeZ(parameters.geometry.lengthZ),
    _stepX(_parameters.bfStep.xRatio * parameters.geometry.lengthX),
    _stepY(_parameters.bfStep.yRatio * parameters.geometry.lengthY),
    _possible(_parameters.parallel.numProcessors[1]*_parameters.parallel.numProcessors[2]==1) { // not working if the domain is split in y and/or z direction

    // for now this post stencil can only be used in sequential simulations or if the domain is only split in x-direction
    // to make it also work in full parallel mode the communication of the uTau values at the wall is necessary

    if (!_possible) return;

    // find the j-index of the first fluid cells on top of the step
    for (int j = 2; j < 2 + _parameters.parallel.localSize[1]; j++) {
      if (_stepY < _parameters.meshsize->getPosY(0,j,0) + 0.5 * _parameters.meshsize->getDy(0,j,0)){
        _jOnStep = j;
        break;
      }
    }

    // find the i-index of the first fluid cells after the step
    _iBehindStep = -1; // indication that the step is not in this domain
    if (_parameters.meshsize->getPosX(1,0,0) < _stepX){
      // only get the index of the cells behind the step if the step is in this domain
      for (int i = 2; i < 2 + _parameters.parallel.localSize[1]; i++) {
        if (_stepX < _parameters.meshsize->getPosX(i,0,0) + 0.5 * _parameters.meshsize->getDx(i,0,0)){
          _iBehindStep = i;
          break;
        }
      }
    }

    const bool is3d = _parameters.geometry.dim == 3;
    if (is3d){
      // sizeX by sizeZ
      _uTau_bottom = (FLOAT*) malloc(_parameters.parallel.localSize[0] * _parameters.parallel.localSize[2] * sizeof(FLOAT));
      _uTau_top    = (FLOAT*) malloc(_parameters.parallel.localSize[0] * _parameters.parallel.localSize[2] * sizeof(FLOAT));
      // sizeX by sizeY
      _uTau_front  = (FLOAT*) malloc(_parameters.parallel.localSize[0] * _parameters.parallel.localSize[1] * sizeof(FLOAT));
      _uTau_back   = (FLOAT*) malloc(_parameters.parallel.localSize[0] * _parameters.parallel.localSize[1] * sizeof(FLOAT));

      // _jOnStep by sizeZ
      _uTau_step   = (FLOAT*) malloc((_jOnStep - 2) * _parameters.parallel.localSize[2] * sizeof(FLOAT));
    }else{
      // sizeX by 1
      _uTau_bottom = (FLOAT*) malloc(_parameters.parallel.localSize[0] * sizeof(FLOAT));
      _uTau_top    = (FLOAT*) malloc(_parameters.parallel.localSize[0] * sizeof(FLOAT));

      // _jOnStep by 1
      _uTau_step   = (FLOAT*) malloc((_jOnStep - 2) * sizeof(FLOAT));
    }
}

TurbulentWallPostStencil::~TurbulentWallPostStencil(){
    if (!_possible) return;

    free(_uTau_bottom);
    free(_uTau_top);
    if (_parameters.geometry.dim == 3){
      free(_uTau_front);
      free(_uTau_back);
    }
    free(_uTau_step);
}

void TurbulentWallPostStencil::preapply ( TurbulentFlowField & turbFlowField ) {
    if (!_possible) return;

    // this preapply method calculates the uTau values at all walls

    Meshsize *ms = _parameters.meshsize;
    FLOAT* velocity = new FLOAT(3);
    int jBottom;
    const int jTop = turbFlowField.getNy() + 1;
    const int kFront = 2, kBack = turbFlowField.getNz() + 1;

    for (int i = 2; i < 2 + turbFlowField.getNx(); i++) {
      if (ms->getPosX(i,0,0) < _stepX){
        jBottom = _jOnStep;
      }else { jBottom = 2; }

      for (int k = 2; k < 2 + turbFlowField.getNz(); k++) {
        turbFlowField.getVelocityCenter(velocity, i, jBottom, k);
        _uTau_bottom[(k-2)*turbFlowField.getNx()+(i-2)] = sqrt(sqrt(pow(velocity[0],2) + pow(velocity[2],2)) / (0.5 * ms->getDy(i,jBottom,k)) / _parameters.flow.Re);

        turbFlowField.getVelocityCenter(velocity, i, jTop, k);
        _uTau_top[(k-2)*turbFlowField.getNx()+(i-2)] = sqrt(sqrt(pow(velocity[0],2) + pow(velocity[2],2)) / (0.5 * ms->getDy(i,jTop,k)) / _parameters.flow.Re);
      }
    }

    for (int i = 2; i < 2 + turbFlowField.getNx(); i++) {
      for (int j = 2; j < 2 + turbFlowField.getNy(); j++) {
        turbFlowField.getVelocityCenter(velocity, i, j, kFront);
        _uTau_front[(j-2)*turbFlowField.getNx()+(i-2)] = sqrt(sqrt(pow(velocity[0],2) + pow(velocity[1],2)) / (0.5 * ms->getDz(i,j,kFront)) / _parameters.flow.Re);

        turbFlowField.getVelocityCenter(velocity, i, j, kBack);
        _uTau_back[(j-2)*turbFlowField.getNx()+(i-2)] = sqrt(sqrt(pow(velocity[0],2) + pow(velocity[1],2)) / (0.5 * ms->getDz(i,j,kBack)) / _parameters.flow.Re);
      }
    }

    if (_iBehindStep > -1){
      for (int j = 2; j < _jOnStep; j++){
        for (int k = 2; k < 2 + turbFlowField.getNz(); k++) {
          turbFlowField.getVelocityCenter(velocity, _iBehindStep, j, k);
          _uTau_step[(j-2)*turbFlowField.getNz()+(k-2)] = sqrt(sqrt(pow(velocity[1],2) + pow(velocity[2],2)) / (0.5 * ms->getDx(_iBehindStep,j,k)) / _parameters.flow.Re);
        }
      }
    }
}

void TurbulentWallPostStencil::apply ( TurbulentFlowField & turbFlowField, int i, int j ) {
  // not yet implemented for 2D
  // TODO: is it necessary for 2D? will we need it?
}

void TurbulentWallPostStencil::apply ( TurbulentFlowField & turbFlowField, int i, int j, int k ) {
    if (!_possible) return;

    Meshsize *ms = _parameters.meshsize;
    FLOAT uTau = 0.0;
    FLOAT yPlus = 0.0;
    FLOAT uPlus[3];

    const int obstacle = turbFlowField.getFlags().getValue(i,j,k);

    // check whether current cell is a fluid cell or not
    if ((obstacle & OBSTACLE_SELF) == 0) {
      turbFlowField.getVelocityCenter(uPlus, i, j, k); // not yet uPlus but normal velocity

      const FLOAT posX = ms->getPosX(i,j,k) + 0.5*ms->getDx(i,j,k);
      const FLOAT posY = ms->getPosY(i,j,k) + 0.5*ms->getDy(i,j,k);
      const FLOAT posZ = ms->getPosZ(i,j,k) + 0.5*ms->getDz(i,j,k);

      // the minimum distance to the wall was already calculated and is assumed to be correct
      // hence, it can be used to find the closest wall
      const FLOAT minDist = turbFlowField.getDistNearestWall().getScalar(i,j,k);

      // if more than one wall is closest, take the maximum uTau value to produce also the maximum yPlus value
      FLOAT uTauWall = -MY_FLOAT_MAX;
      if ((_sizeY - posY) <= minDist){
        // distance is minimal to the top
        uTauWall = std::max(uTauWall, _uTau_top[(k-2)*turbFlowField.getNx()+(i-2)]);
      }
      if ((posX < _stepX && (posY - _stepY) <= minDist) || (_stepX < posX && posY <= minDist)){
        // distance is minimal to the bottom
        uTauWall = std::max(uTauWall, _uTau_bottom[(k-2)*turbFlowField.getNx()+(i-2)]);
      }
      if ((_sizeZ - posZ) <= minDist){
        // distance is minimal to the back
        uTauWall = std::max(uTauWall, _uTau_back[(j-2)*turbFlowField.getNx()+(i-2)]);
      }
      if (posZ <= minDist){
        // distance is minimal to the front
        uTauWall = std::max(uTauWall, _uTau_front[(j-2)*turbFlowField.getNx()+(i-2)]);
      }
      if (_stepX < posX && posY < _stepY && (posX - _stepX) <= minDist){
        // distance is minimal to the back face of the step
        uTauWall = std::max(uTauWall, _uTau_step[(j-2)*turbFlowField.getNz()+(k-2)]);
      }

      if ((j == 2 && posX > _stepX) || (j == _jOnStep && posX < _stepX) || (i == _iBehindStep && posY < _stepY)
          || (k == 2) || (j == turbFlowField.getNy() + 1) || (k == turbFlowField.getNz() + 1)){

        uTau = uTauWall;
      }

      if (posX < _stepX && posY < _stepY)
        std::cout << "??? " << i << " " << j << " " << k << " " << obstacle << std::endl;

      yPlus = 0.5 * minDist * uTauWall * _parameters.flow.Re;

      if (uTauWall > 0.0){
        uPlus[0] = uPlus[0] / uTauWall;
        uPlus[1] = uPlus[1] / uTauWall;
        uPlus[2] = uPlus[2] / uTauWall;
      }else{
        uPlus[0] = 0.0;
        uPlus[1] = 0.0;
        uPlus[2] = 0.0;
      }
    } else {
      uPlus[0] = 0.0;
      uPlus[1] = 0.0;
      uPlus[2] = 0.0;
    }

    _ssUTau << uTau << std::endl;
    _ssYPlus << yPlus << std::endl;
    _ssUPlus << uPlus[0] << " " << uPlus[1] << " " << uPlus[2] << std::endl;
}

void TurbulentWallPostStencil::writeVtk ( std::ostream & stream ) {
    if (!_possible) return;

    // write wall shear velocity data
    stream << "\nSCALARS uTau float 1" << std::endl;
    stream << "LOOKUP_TABLE default" << std::endl;
    stream << _ssUTau.str();

    // write non-dimensional wall distance data
    stream << "\nSCALARS yPlus float 1" << std::endl;
    stream << "LOOKUP_TABLE default" << std::endl;
    stream << _ssYPlus.str();

    // write non dimensional velocity
    stream << "\nVECTORS uPlus float" << std::endl;
    stream << _ssUPlus.str();

    // clear the stringstreams
    _ssUTau.str("");
    _ssUPlus.str("");
    _ssYPlus.str("");
}
