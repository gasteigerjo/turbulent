#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "VTKStencil.h"


VTKStencil::VTKStencil ( const Parameters & parameters) : FieldStencil<FlowField> (parameters), _prefix (parameters.vtk.prefix) {

    // Get number of cells (i.e. array length)
    int numCells = 1;
    for(int i = 0; i < _parameters.geometry.dim; i++) {
        numCells *= _parameters.parallel.localSize[i];
    }

    // Set up arrays for temporarily saving the pressure and velocity
    _pressures = new FLOAT[numCells];
    _velocities = new FLOAT[3*numCells];
}


void VTKStencil::apply ( FlowField & flowField, int i, int j ){
    if( i > 1 && j > 1 ) { // only for fluid cells

        //temporarily save pressure and interpolated velocity in an array
        int ind = (j-2) * _parameters.parallel.localSize[0] + (i-2);
        flowField.getPressureAndVelocity(_pressures[ind], _velocities + 3*ind, i, j);
    }
}


void VTKStencil::apply ( FlowField & flowField, int i, int j, int k ){
    if( i > 1 && j > 1 && k > 1 ) { // only for fluid cells

        //temporarily save pressure and interpolated velocity in an array
        int ind = ((k-2) * _parameters.parallel.localSize[1] + (j-2)) * _parameters.parallel.localSize[0] + (i-2);
        flowField.getPressureAndVelocity(_pressures[ind], _velocities + 3*ind, i, j, k);
    }
}

void VTKStencil::write ( FlowField & flowField, int timeStep ){

    // set up variables for the number of cells
    int numCells = 1;
    int numVertices = 1;
    int cellsMin[3];
    int cellsMax[3];
    for(int i = 0; i < 3; i++) {
        cellsMin[i] = _parameters.parallel.firstCorner[i];
        cellsMax[i] = cellsMin[i] + _parameters.parallel.localSize[i];
        if (i < _parameters.geometry.dim) {
            numCells *= _parameters.parallel.localSize[i];
            numVertices *= _parameters.parallel.localSize[i] + 1;
        }
    }

    // open stream
    std::ofstream vtkFile;
    std::stringstream filename;
    filename << _prefix << "_" << std::setfill('0') << std::setw(6) << timeStep << ".vtk";
    vtkFile.open(filename.str().c_str());
    if (vtkFile.fail()) {
        std::cout << "Failed to open file for VTK output: " << filename.str() << std::endl;
    } else {
        std::cout << "Writing VTK output: " << filename.str() << std::endl;
    }

    // write header
    vtkFile << "# vtk DataFile Version 2.0" << std::endl << "I need something to put here" << std::endl << "ASCII" << std::endl << std::endl;

    // write vertex coordinates
    vtkFile << "DATASET STRUCTURED_GRID" << std::endl << "DIMENSIONS "
        << _parameters.parallel.localSize[0] + 1 << " "
        << _parameters.parallel.localSize[1] + 1 << " "
        << _parameters.parallel.localSize[2] + 1 << std::endl
        << "POINTS " << numVertices << " float" << std::endl;

    vtkFile << std::fixed << std::setprecision(6);
    for(int k = cellsMin[2]; k <= cellsMax[2]; k++) {
        for(int j = cellsMin[1]; j <= cellsMax[1]; j++) {
            for(int i = cellsMin[0]; i <= cellsMax[0]; i++) {
                vtkFile << _parameters.meshsize->getPosX(i+2, j+2, k+2) << " "
                    << _parameters.meshsize->getPosY(i+2, j+2, k+2) << " "
                    << _parameters.meshsize->getPosZ(i+2, j+2, k+2) << std::endl;
            }
        }
    }
    vtkFile << std::endl;

    // write pressure
    int ind;
    vtkFile << "CELL_DATA " << numCells << std::endl << "SCALARS pressure float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    vtkFile << std::scientific;
    for(int k = cellsMin[2]; k < cellsMax[2] + (3-_parameters.geometry.dim); k++) {
        for(int j = cellsMin[1]; j < cellsMax[1]; j++) {
            for(int i = cellsMin[0]; i < cellsMax[0]; i++) {
                ind = (k * _parameters.parallel.localSize[1] + j) * _parameters.parallel.localSize[0] + i;
                vtkFile << _pressures[ind] << std::endl;
            }
        }
    }
    vtkFile << std::endl;

    // write velocity
    vtkFile << "VECTORS velocity float" << std::endl;
    for(int k = cellsMin[2]; k < cellsMax[2] + (3-_parameters.geometry.dim); k++) {
        for(int j = cellsMin[1]; j < cellsMax[1]; j++) {
            for(int i = cellsMin[0]; i < cellsMax[0]; i++) {
                ind = (k * _parameters.parallel.localSize[1] + j) * _parameters.parallel.localSize[0] + i;
                vtkFile << _velocities[3*ind] << " " << _velocities[3*ind + 1] << " " << _velocities[3*ind + 2] << std::endl;
            }
        }
    }
    vtkFile << std::endl;

    vtkFile.close();
}
