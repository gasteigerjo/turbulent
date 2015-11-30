#include "VTKStencil.h"
#include "../Iterators.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "../TurbulentFlowField.h"

VTKStencil::VTKStencil ( const Parameters & parameters ) :
  FieldStencil<FlowField> ( parameters ),
  _turbulent(parameters.simulation.type=="turbulence")
  {}


void VTKStencil::apply ( FlowField & flowField, int i, int j ) {
    Meshsize *ms = _parameters.meshsize;
    _ssPoints << ms->getPosX(i+1, j+1) << " " << ms->getPosY(i+1, j+1) << " 0" << std::endl;

    // skip ghost cells
    if (i < 2 || j < 2) return;

    FLOAT pressure, turbVisc;
    FLOAT* velocity = new FLOAT(2);
    const int obstacle = flowField.getFlags().getValue(i, j);

    // check whether current cell is a fluid cell or not
    if ((obstacle & OBSTACLE_SELF) == 0) {
      flowField.getPressureAndVelocity(pressure, velocity,  i, j);
      if (_turbulent){
        TurbulentFlowField* turbFlowField = static_cast<TurbulentFlowField*>(&flowField);
        turbVisc = turbFlowField->getTurbViscosity().getScalar(i, j);
        // std::cout << turbVisc;
      }
    } else {
      pressure = 0.0;
      velocity[0] = 0.0;
      velocity[1] = 0.0;
      turbVisc = 0.0;
    }

    _ssPressure << pressure << std::endl;
    _ssVelocity << velocity[0] << " " << velocity[1] << " 0" << std::endl;
    if (_turbulent) _ssTurbViscosity << turbVisc << std::endl;
}


void VTKStencil::apply ( FlowField & flowField, int i, int j, int k ) {
    Meshsize *ms = _parameters.meshsize;
    _ssPoints << ms->getPosX(i+1, j+1, k+1) << " "
              << ms->getPosY(i+1, j+1, k+1) << " "
              << ms->getPosZ(i+1, j+1, k+1) << std::endl;

    // skip ghost cells
    if (i < 2 || j < 2 || k < 2) return;

    FLOAT pressure, turbVisc;
    FLOAT* velocity = new FLOAT(3);
    const int obstacle = flowField.getFlags().getValue(i, j, k);

    // check whether current cell is a fluid cell or not
    if ((obstacle & OBSTACLE_SELF) == 0) {
      flowField.getPressureAndVelocity(pressure, velocity,  i, j, k);
      if (_turbulent){
        TurbulentFlowField* turbFlowField = static_cast<TurbulentFlowField*>(&flowField);
        turbVisc = turbFlowField->getTurbViscosity().getScalar(i, j, k);
      }
    } else {
      pressure = 0.0;
      velocity[0] = 0.0;
      velocity[1] = 0.0;
      velocity[2] = 0.0;
      turbVisc = 0.0;
    }

    _ssPressure << pressure << std::endl;
    _ssVelocity << velocity[0] << " " << velocity[1] << " " << velocity[2] << std::endl;
    if (_turbulent) _ssTurbViscosity << turbVisc << std::endl;
}

void VTKStencil::write ( FlowField & flowField, int timeStep ) {
    const int nx = flowField.getNx();
    const int ny = flowField.getNy();
    const int nz = _parameters.geometry.dim == 2 ? 0 : flowField.getNz();
    const int nPoints = (nx+1)*(ny+1)*(nz+1);
    const int nCells =_parameters.geometry.dim == 2 ? nx*ny :  nx*ny*nz;

    // construct the file name and open the corresponding file
    std::stringstream fileName;
    fileName << _parameters.vtk.prefix << "_" << std::setfill('0') << std::setw(6) << timeStep << ".vtk";
    std::ofstream file;
    file.open(fileName.str().c_str());
    fileName.clear();

    // write VTK header
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "I need something to put here" << std::endl;
    file << "ASCII\n" << std::endl;

    // write grid
    file << "DATASET STRUCTURED_GRID" << std::endl;
    file << "DIMENSIONS " << nx+1 << " " << ny+1 << " " << nz+1 << std::endl;
    file << "POINTS " << nPoints << " float" << std::endl;
    file << _ssPoints.str();

    file << "\nCELL_DATA " << nCells << std::endl;

    // write pressure data
    file << "SCALARS pressure float 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    file << _ssPressure.str();

    // write velocity data
    file << "\nVECTORS velocity float" << std::endl;
    file << _ssVelocity.str();

    if (_turbulent){
      // write pressure data
      file << "\nSCALARS turbViscosity float 1" << std::endl;
      file << "LOOKUP_TABLE default" << std::endl;
      file << _ssTurbViscosity.str();
    }


    // clear the stringstreams
    _ssPoints.str("");
    _ssPressure.str("");
    _ssVelocity.str("");
    _ssTurbViscosity.str("");

    // finally, close the file
    file.close();
}
