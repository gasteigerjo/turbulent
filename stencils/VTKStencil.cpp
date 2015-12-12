#include "VTKStencil.h"
#include "../Iterators.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "../TurbulentFlowField.h"

VTKStencil::VTKStencil ( const Parameters & parameters ) :
  FieldStencil<FlowField> ( parameters ),
  _turbulent(parameters.simulation.type=="turbulence"),
  _includeGhostCells(false)
  {}

VTKStencil::VTKStencil ( const Parameters & parameters, bool includeGhostCells ) :
  FieldStencil<FlowField> ( parameters ),
  _turbulent(parameters.simulation.type=="turbulence"),
  _includeGhostCells(includeGhostCells)
  {}


void VTKStencil::apply ( FlowField & flowField, int i, int j ) {
    Meshsize *ms = _parameters.meshsize;
    // Sizes per x and y direction, including the three ghost layers of each
    const int cellsX = flowField.getCellsX();
    const int cellsY = flowField.getCellsY();

    // Ghost cells only
    // TODO: confirm that correct vtk files are produced (correct order of coordinates)
    if (_includeGhostCells && i==0){
      if (j==0){
        for (int k = 0; k < cellsX+1; k++) {
          _ssPoints << ms->getPosX(k, 0) << " " << ms->getPosY(k, 0) << " 0" << std::endl;
        }
      }
      _ssPoints << ms->getPosX(0, j+1) << " " << ms->getPosY(0, j+1) << " 0" << std::endl;
    }
    // NOTE: we need also the first layer to construct CELL data (confirm the implementation)
    if (_includeGhostCells || (i > 1 && j > 1 && i < cellsX && j < cellsY)) {
      _ssPoints << ms->getPosX(i, j) << " " << ms->getPosY(i, j) << " 0" << std::endl;
    }

    // skip assigning property values to ghost cells, if not asked to include them.
    if (!_includeGhostCells && (i < 2 || j < 2 || i > cellsX-2 || j > cellsY-2)) return;

    FLOAT pressure, turbVisc;
    // FLOAT distWall;
    FLOAT* velocity = new FLOAT(2);
    const int obstacle = flowField.getFlags().getValue(i, j);

    // check whether current cell is a fluid cell or not
    if ((obstacle & OBSTACLE_SELF) == 0) {
      flowField.getPressureAndVelocity(pressure, velocity,  i, j);
      if (_turbulent){
        TurbulentFlowField* turbFlowField = static_cast<TurbulentFlowField*>(&flowField);
        turbVisc = turbFlowField->getTurbViscosity().getScalar(i, j);
        // distWall = turbFlowField->getDistNearestWall().getScalar(i, j);
      }
    } else {
      pressure = 0.0;
      velocity[0] = 0.0;
      velocity[1] = 0.0;
      turbVisc = 0.0;
      // distWall = 0.0;
    }

    _ssPressure << pressure << std::endl;
    _ssVelocity << velocity[0] << " " << velocity[1] << " 0" << std::endl;
    _ssFlags << obstacle << std::endl;
    if (_turbulent){
      _ssTurbViscosity << turbVisc << std::endl;
      // _ssDistWall << distWall << std::endl;
    }
}

void VTKStencil::apply ( FlowField & flowField, int i, int j, int k ) {
    Meshsize *ms = _parameters.meshsize;
    // Sizes per x, y and z direction, including the three ghost layers of each
    const int cellsX = flowField.getCellsX();
    const int cellsY = flowField.getCellsY();
    const int cellsZ = flowField.getCellsZ();

    // Ghost cells only
    // TODO: implement the ghost layer inclusion in 3D.

    if (_includeGhostCells || (i > 1 && j > 1 && k > 1 && i < cellsX && j < cellsY && k < cellsZ)) {
      _ssPoints << ms->getPosX(i, j, k) << " " << ms->getPosY(i, j, k) << " " << ms->getPosZ(i, j, k) << std::endl;
    }

    // skip assigning property values to ghost cells, if not asked to include them.
    if (!_includeGhostCells && (i < 2 || j < 2 || k < 2 || i > cellsX-2 || j > cellsY-2 || k > cellsZ-2 )) return;

    FLOAT pressure, turbVisc;
    // FLOAT distWall;
    FLOAT* velocity = new FLOAT(3);
    const int obstacle = flowField.getFlags().getValue(i, j, k);
    // FLOAT* fgh = flowField.getFGH().getVector(i,j,k);

    // check whether current cell is a fluid cell or not
    if ((obstacle & OBSTACLE_SELF) == 0) {
      flowField.getPressureAndVelocity(pressure, velocity,  i, j, k);
      if (_turbulent){
        TurbulentFlowField* turbFlowField = static_cast<TurbulentFlowField*>(&flowField);
        turbVisc = turbFlowField->getTurbViscosity().getScalar(i, j, k);
        // distWall = turbFlowField->getDistNearestWall().getScalar(i, j, k);
      }
    } else {
      pressure = 0.0;
      velocity[0] = 0.0;
      velocity[1] = 0.0;
      velocity[2] = 0.0;
      turbVisc = 0.0;
      // distWall = 0.0;
    }

    _ssPressure << pressure << std::endl;
    _ssVelocity << velocity[0] << " " << velocity[1] << " " << velocity[2] << std::endl;
    _ssFlags << obstacle << std::endl;
    if (_turbulent){
      _ssTurbViscosity << turbVisc << std::endl;
      // _ssDistWall << distWall << std::endl;
    }
    // _ssFGH << fgh[0] << " " << fgh[1] << " " << fgh[2] << std::endl;
}

void VTKStencil::write ( FlowField & flowField, int timeStep ) {
    int nx = flowField.getNx();
    int ny = flowField.getNy();
    int nz = _parameters.geometry.dim == 2 ? 0 : flowField.getNz();
    if (_includeGhostCells){
      nx = flowField.getCellsX();
      ny = flowField.getCellsY();
      nz = _parameters.geometry.dim == 2 ? 0 : flowField.getCellsZ();
    }
    const int nPoints = (nx+1)*(ny+1)*(nz+1);
    const int nCells =_parameters.geometry.dim == 2 ? nx*ny :  nx*ny*nz;

    // construct the file name and open the corresponding file
    std::stringstream fileName;
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
    file << _ssPoints.str();

    file << "\nCELL_DATA " << nCells << std::endl;

    // write pressure data
    file << "\nSCALARS pressure float 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    file << _ssPressure.str();

    // write velocity data
    file << "\nVECTORS velocity float" << std::endl;
    file << _ssVelocity.str();

    // // write FGH data
    // file << "\nVECTORS fgh float" << std::endl;
    // file << _ssFGH.str();

    // // write flags data
    // file << "\nSCALARS flags integer 1" << std::endl;
    // file << "LOOKUP_TABLE default" << std::endl;
    // file << _ssFlags.str();

    if (_turbulent){
      // write turbulent viscosity data
      file << "\nSCALARS turbViscosity float 1" << std::endl;
      file << "LOOKUP_TABLE default" << std::endl;
      file << _ssTurbViscosity.str();

      // // write wall distance data
      // file << "\nSCALARS distWall float 1" << std::endl;
      // file << "LOOKUP_TABLE default" << std::endl;
      // file << _ssDistWall.str();

      // clear the turbulent specific stringstreams
      _ssTurbViscosity.str("");
      // _ssDistWall.str("");
    }


    // clear the stringstreams
    _ssPoints.str("");
    _ssPressure.str("");
    _ssVelocity.str("");
    // _ssFGH.str("");
    // _ssFlags.str("");

    // finally, close the file
    file.close();
}
