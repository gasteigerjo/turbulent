#include "Meshsize.h"

#include "Parameters.h"

UniformMeshsize::UniformMeshsize(
  const Parameters &parameters
): Meshsize(),
_dx(parameters.geometry.lengthX/parameters.geometry.sizeX),
_dy(parameters.geometry.lengthY/parameters.geometry.sizeY),
_dz(parameters.geometry.dim==3 ? parameters.geometry.lengthZ/parameters.geometry.sizeZ : 0.0),
_firstCornerX(parameters.parallel.firstCorner[0]),
_firstCornerY(parameters.parallel.firstCorner[1]),
_firstCornerZ(parameters.geometry.dim==3 ? parameters.parallel.firstCorner[2] : 0)
{
  if (_dx<= 0.0){ handleError(1,"_dx<=0.0!"); }
  if (_dy<= 0.0){ handleError(1,"_dy<=0.0!"); }
  if (parameters.geometry.dim==3){
    if (_dz<=0.0){ handleError(1,"_dz<=0.0!"); }
  }
}

UniformMeshsize::~UniformMeshsize(){}



TanhMeshStretching::TanhMeshStretching(
  const Parameters & parameters,bool stretchX, bool stretchY, bool stretchZ
): Meshsize(), _parameters(parameters), _uniformMeshsize(parameters),
   _lengthX(parameters.geometry.lengthX), _lengthY(parameters.geometry.lengthY),
   _lengthZ(parameters.geometry.dim==3 ? parameters.geometry.lengthZ : 0.0),
   _sizeX(parameters.geometry.sizeX), _sizeY(parameters.geometry.sizeY),
   _sizeZ(parameters.geometry.dim==3 ? parameters.geometry.sizeZ : 1),
   _firstCornerX(parameters.parallel.firstCorner[0]), _firstCornerY(parameters.parallel.firstCorner[1]),
   _firstCornerZ(parameters.geometry.dim==3 ? parameters.parallel.firstCorner[2] : 0),
   _stretchX(stretchX), _stretchY(stretchY), _stretchZ(stretchZ),
   _dxMin(stretchX ? 0.5*parameters.geometry.lengthX*(1.0 + tanh(parameters.geometry.deltaSX*(2.0/parameters.geometry.sizeX-1.0))/tanh(parameters.geometry.deltaSX)) : _uniformMeshsize.getDx(0,0,0)),
   _dyMin(stretchY ? 0.5*parameters.geometry.lengthY*(1.0 + tanh(parameters.geometry.deltaSY*(2.0/parameters.geometry.sizeY-1.0))/tanh(parameters.geometry.deltaSY)) : _uniformMeshsize.getDy(0,0,0)),
   _dzMin(stretchZ ? 0.5*parameters.geometry.lengthZ*(1.0 + tanh(parameters.geometry.deltaSZ*(2.0/parameters.geometry.sizeZ-1.0))/tanh(parameters.geometry.deltaSZ)) : _uniformMeshsize.getDz(0,0,0))
{ }

TanhMeshStretching::~TanhMeshStretching(){
  free(_coordinatesX);
  free(_coordinatesY);
  free(_coordinatesZ);
}

void TanhMeshStretching::precomputeCoordinates(){
  // + 4 is needed to also get the height of the last cell (+ 3 due to the three ghost layers)
  if (_stretchX){
    _coordinatesX = (FLOAT*) malloc((_parameters.parallel.localSize[0] + 4) * sizeof(FLOAT));
    for (int i = 0; i < _parameters.parallel.localSize[0] + 4; i++) {
      _coordinatesX[i] = computeCoordinateX(i);
    }
  }
  if (_stretchY){
    _coordinatesY = (FLOAT*) malloc((_parameters.parallel.localSize[1] + 4) * sizeof(FLOAT));
    for (int i = 0; i < _parameters.parallel.localSize[1] + 4; i++) {
      _coordinatesY[i] = computeCoordinateY(i);
    }
  }
  if (_stretchZ){
    _coordinatesZ = (FLOAT*) malloc((_parameters.parallel.localSize[2] + 4) * sizeof(FLOAT));
    for (int i = 0; i < _parameters.parallel.localSize[2] + 4; i++) {
      _coordinatesZ[i] = computeCoordinateZ(i);
    }
  }
}

FLOAT TanhMeshStretching::computeCoordinate(int i, int firstCorner, int size, FLOAT length, FLOAT deltaS, FLOAT dxMin) const {
  const int index = i-2+firstCorner;
  // equidistant mesh on lower/left part
  if (index < 0){
    return dxMin*index;
  // equidistant mesh on upper/right part
  } else if (index > size-1){
    return length+dxMin*(index-size);
  } else {
    // stretched mesh on lower half of channel -> we check if we are in lower 50% and then use stretching for 2.0*p
    FLOAT p = ((FLOAT) index)/size;
    if (p<0.5){
      return 0.5*length*(1.0 + tanh(deltaS*(2.0*p-1.0))/tanh(deltaS));
    // stretched mesh on upper half of channel -> we mirror the stretching
    } else {
      p = ((FLOAT) size-index)/size;
      return length-0.5*length*(1.0 + tanh(deltaS*(2.0*p-1.0))/tanh(deltaS));
    }
  }
}

FLOAT TanhMeshStretching::computeCoordinateX(int i) const {
  return computeCoordinate(i,_firstCornerX,_sizeX,_lengthX,_parameters.geometry.deltaSX,_dxMin);
}

FLOAT TanhMeshStretching::computeCoordinateY(int i) const {
  return computeCoordinate(i,_firstCornerY,_sizeY,_lengthY,_parameters.geometry.deltaSY,_dyMin);
}

FLOAT TanhMeshStretching::computeCoordinateZ(int i) const {
  return computeCoordinate(i,_firstCornerZ,_sizeZ,_lengthZ,_parameters.geometry.deltaSZ,_dzMin);
}



BfsMeshStretching::BfsMeshStretching(
  const Parameters & parameters
): TanhMeshStretching(parameters,true,true,true), _stepX(parameters.bfStep.xRatio * parameters.geometry.lengthX)
{ }

BfsMeshStretching::~BfsMeshStretching(){}

void BfsMeshStretching::precomputeCoordinates(){
  FLOAT deltaS = _parameters.geometry.deltaSX;
  _sizeXBeforeStep = _parameters.bfStep.xRatio * (_sizeX - 4) + 2;
  _sizeXAfterStep = _sizeX - _sizeXBeforeStep;
  _dxMinBefore = _stepX/(_sizeXBeforeStep - 2)*(1.0 + tanh(deltaS*(1.0/3.0-1.0))/tanh(deltaS));
  _dxMinAfter = (_lengthX - _stepX)/(_sizeXAfterStep - 2)*(1.0 + tanh(deltaS*(1.0/3.0-1.0))/tanh(deltaS));

  deltaS = _parameters.geometry.deltaSY;
  _dyMinBelow = 1.0;
  _dyMinAbove = 0.0;
  _sizeYBelowStep = _parameters.bfStep.yRatio * _sizeY - 1;
  _sizeYAboveStep = _sizeY - _sizeYBelowStep;
  while(_dyMinAbove < _dyMinBelow) {
    _sizeYBelowStep++;
    _sizeYAboveStep--;
    _dyMinBelow = 0.5*_parameters.bfStep.yRatio*_lengthY*(1.0 + tanh(deltaS*(2.0/_sizeYBelowStep-1.0))/tanh(deltaS));
    _dyMinAbove = 0.5*(1.0-_parameters.bfStep.yRatio)*_lengthY*(1.0 + tanh(deltaS*(2.0/_sizeYAboveStep-1.0))/tanh(deltaS));
  }

  TanhMeshStretching::precomputeCoordinates();
}

FLOAT BfsMeshStretching::computeCoordinateX(int i) const {
  const FLOAT deltaS = _parameters.geometry.deltaSX;
  const int index = i - 2 + _firstCornerX;
  if (index < _sizeXBeforeStep){
    const FLOAT dx = _stepX / (_sizeXBeforeStep - 2);
    if (index <= _sizeXBeforeStep - 3){
      return index * dx;
    }else{
      return _stepX - dx*(1.0 + tanh(deltaS*((_sizeXBeforeStep-index)/3.0-1.0))/tanh(deltaS));
    }
  }else if (index == _sizeXBeforeStep){
    return _stepX;
  }else{
    const FLOAT dx = (_lengthX - _stepX) / (_sizeXAfterStep - 2);
    if (index < _sizeXBeforeStep + 3){
      return _stepX +  dx*(1.0 + tanh(deltaS*((index-_sizeXBeforeStep)/3.0-1.0))/tanh(deltaS));
    }else{
      return _stepX + (index - _sizeXBeforeStep - 2) * dx;
    }
  }
}

FLOAT BfsMeshStretching::computeCoordinateY(int i) const {
  const int index = i - 2 + _firstCornerY;
  if (index < _sizeYBelowStep){
    return computeCoordinate(i,_firstCornerY,_sizeYBelowStep,_parameters.bfStep.yRatio*_lengthY,_parameters.geometry.deltaSY,_dyMinBelow);
  }else{
    return _parameters.bfStep.yRatio*_lengthY + computeCoordinate(i-_sizeYBelowStep,_firstCornerY,_sizeYAboveStep,(1.0-_parameters.bfStep.yRatio)*_lengthY,_parameters.geometry.deltaSY,_dyMinAbove);
  }
}
