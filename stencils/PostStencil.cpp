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
