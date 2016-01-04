#include "VtkOutput.h"
#include "Iterators.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TurbulentFlowField.h"
#include "stencils/PostStencil.h"
#include "stencils/TurbulentPostStencil.h"

VtkOutput::VtkOutput ( FlowField & flowField, const Parameters & parameters ) :
  FieldStencil<FlowField>(parameters),
  _turbulent(parameters.simulation.type=="turbulence"),
  _flowField(flowField),
  _postIterator(flowField,parameters,*this,1,0)
{
    _nPostStencils = 1;
    _postStencils = new PostStencil<FlowField>*[_nPostStencils];
    _postStencils[0] = new BasicPostStencil(_parameters);

    if (_turbulent){
      _turbFlowField = (TurbulentFlowField*) &flowField;

      _nTurbPostStencils = 2;
      _turbPostStencils = new PostStencil<TurbulentFlowField>*[_nTurbPostStencils];
      _turbPostStencils[0] = new TurbulentPostStencil(_parameters);
      _turbPostStencils[1] = new TurbulentWallPostStencil(_parameters);
    }
}

VtkOutput::~VtkOutput (){
  free(_postStencils);
  free(_turbPostStencils);
}


void VtkOutput::apply ( FlowField & flowField, int i, int j ) {
    // apply all post stencils
    // (the implementation of this class as a FieldStencil allows the use of
    // a single iterator and this method to easily iterate all post stencils)

    for (int post = 0; post < _nPostStencils; post++) {
      _postStencils[post]->apply(_flowField,i,j);
    }

    if (_turbulent){
      for (int post = 0; post < _nTurbPostStencils; post++) {
        _turbPostStencils[post]->apply(*_turbFlowField,i,j);
      }
    }
}


void VtkOutput::apply ( FlowField & flowField, int i, int j, int k ) {
    // apply all post stencils
    // (the implementation of this class as a FieldStencil allows the use of
    // a single iterator and this method to easily iterate all post stencils)

    for (int post = 0; post < _nPostStencils; post++) {
      _postStencils[post]->apply(_flowField,i,j,k);
    }

    if (_turbulent){
      for (int post = 0; post < _nTurbPostStencils; post++) {
        _turbPostStencils[post]->apply(*_turbFlowField,i,j,k);
      }
    }
}

void VtkOutput::write ( int timeStep ) {
    // preapply the post stencils
    for (int post = 0; post < _nPostStencils; post++) {
      _postStencils[post]->preapply(_flowField);
    }
    // also preapply the turbulent post stencils if existent
    if (_turbulent){
      for (int post = 0; post < _nTurbPostStencils; post++) {
        _turbPostStencils[post]->preapply(*_turbFlowField);
      }
    }

    // iterate all post stencils by iterating "this"
    _postIterator.iterate();

    Meshsize *ms = _parameters.meshsize;
    int nx = _flowField.getNx();
    int ny = _flowField.getNy();
    int nz = _parameters.geometry.dim == 2 ? 0 : _flowField.getNz();
    const int nPoints = (nx+1)*(ny+1)*(nz+1);
    const int nCells = _parameters.geometry.dim == 2 ? nx*ny :  nx*ny*nz;

    // construct the file name and open the corresponding file
    std::ostringstream fileName;
    fileName << _parameters.vtk.prefix
             << "_" << std::setfill('0') << std::setw(4) << _parameters.parallel.rank
             << "." << std::setfill('0') << std::setw(6) << timeStep
             << ".vtk";
    std::ofstream file;
    file.open(fileName.str().c_str());
    fileName.clear();

    // write VTK header
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "NS-EOF output for rank = " << _parameters.parallel.rank << " and timestep = " << timeStep << std::endl;
    file << "ASCII\n" << std::endl;


    // write grid
    file << "DATASET STRUCTURED_GRID" << std::endl;
    file << "DIMENSIONS " << nx+1 << " " << ny+1 << " " << nz+1 << std::endl;
    file << "POINTS " << nPoints << " float" << std::endl;

    if (_parameters.geometry.dim == 2){
        for (int j = 2; j < ny + 3; j++){
            for (int i = 2; i < nx + 3; i++){
                file << ms->getPosX(i, j) << " " << ms->getPosY(i, j) << " 0" << std::endl;
            }
        }
    }else if (_parameters.geometry.dim == 3){
        for (int k = 2; k < nz + 3; k++){
            for (int j = 2; j < ny + 3; j++){
                for (int i = 2; i < nx + 3; i++){
                    file << ms->getPosX(i, j, k) << " " << ms->getPosY(i, j, k) << " " << ms->getPosZ(i, j, k) << std::endl;
                }
            }
        }
    }

    file << "\nCELL_DATA " << nCells << std::endl;

    // write the data from the post stencils
    for (int post = 0; post < _nPostStencils; post++) {
      _postStencils[post]->writeVtk(file);
    }

    // also write the data from the turbulent post stencils
    if (_turbulent){
      for (int post = 0; post < _nTurbPostStencils; post++) {
        _turbPostStencils[post]->writeVtk(file);
      }
    }

    // finally, close the file
    file.close();
}
