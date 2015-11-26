#ifndef _TURBULENT_SIMULATION_H_
#define _TURBULENT_SIMULATION_H_

#include "Simulation.h"
#include "TurbulentFlowField.h"
#include "Iterators.h"
#include "stencils/FGHTurbStencil.h"


class TurbulentSimulation : public Simulation {
  protected:
    TurbulentFlowField &_turbFlowField;

    FGHTurbStencil _fghTurbStencil;
    FieldIterator<TurbulentFlowField> _fghTurbIterator;

  public:
    TurbulentSimulation(Parameters &parameters, TurbulentFlowField &turbFlowField):
      Simulation(parameters, turbFlowField),
      _turbFlowField(turbFlowField),
      _fghTurbStencil(parameters),
      _fghTurbIterator(_turbFlowField,parameters,_fghTurbStencil)
    {}
};

#endif // _TURBULENT_SIMULATION_H_
