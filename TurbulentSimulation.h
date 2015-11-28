#ifndef _TURBULENT_SIMULATION_H_
#define _TURBULENT_SIMULATION_H_

#include "Simulation.h"
#include "TurbulentFlowField.h"
#include "Iterators.h"
#include "stencils/FGHTurbStencil.h"
#include "stencils/TurbViscosityStencil.h"
#include "stencils/DistNearestWallStencil.h"


class TurbulentSimulation : public Simulation {
  protected:
    TurbulentFlowField &_turbFlowField;

    FGHTurbStencil _fghTurbStencil;
    FieldIterator<TurbulentFlowField> _fghTurbIterator;

    TurbViscosityStencil _turbViscStencil;
    FieldIterator<TurbulentFlowField> _turbViscIterator;

  public:
    TurbulentSimulation(Parameters &parameters, TurbulentFlowField &turbFlowField):
      Simulation(parameters, turbFlowField),
      _turbFlowField(turbFlowField),
      _fghTurbStencil(parameters),
      _fghTurbIterator(turbFlowField,parameters,_fghTurbStencil),
      _turbViscStencil(parameters),
      _turbViscIterator(turbFlowField,parameters,_turbViscStencil)
    {}

    void initializeFlowField() {
      // first call initialization of base class
      Simulation::initializeFlowField();

      // calculate the distance to the nearest wall for each cell
      DistNearestWallStencil distStencil(_parameters);
      FieldIterator<TurbulentFlowField> iterator(_turbFlowField,_parameters,distStencil);
      iterator.iterate();
    }

    void solveTimestep(){
      // TODO timestepping for turbulent simulation
    }

  protected:
    void setTimeStep(){
      // TODO find implementation, p.8 of work sheet
    }
};

#endif // _TURBULENT_SIMULATION_H_
