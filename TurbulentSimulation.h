#ifndef _TURBULENT_SIMULATION_H_
#define _TURBULENT_SIMULATION_H_

#include "Simulation.h"
#include "TurbulentFlowField.h"
#include "Iterators.h"
#include "stencils/FGHTurbStencil.h"
#include "stencils/TurbViscosityStencil.h"
#include "stencils/DistNearestWallStencil.h"
#include "stencils/MinDtStencil.h"
#include "stencils/TurbViscosityBoundaryStencil.h"


class TurbulentSimulation : public Simulation {
  protected:
    TurbulentFlowField &_turbFlowField;

    FGHTurbStencil _fghTurbStencil;
    FieldIterator<TurbulentFlowField> _fghTurbIterator;

    TurbViscosityStencil _turbViscStencil;
    FieldIterator<TurbulentFlowField> _turbViscIterator;

    MinDtStencil _minDtStencil;
    FieldIterator<TurbulentFlowField> _minDtIterator;

    GlobalBoundaryIterator<TurbulentFlowField> _wallTurbVisIterator;

  public:
    TurbulentSimulation(Parameters &parameters, TurbulentFlowField &turbFlowField):
      Simulation(parameters, turbFlowField),
      _turbFlowField(turbFlowField),
      _fghTurbStencil(parameters),
      _fghTurbIterator(turbFlowField,parameters,_fghTurbStencil),
      _turbViscStencil(parameters),
      _turbViscIterator(turbFlowField,parameters,_turbViscStencil),
      _minDtStencil(parameters),
      _minDtIterator(turbFlowField,parameters,_minDtStencil),
      _wallTurbVisIterator(createGlobalBoundaryTurbViscIterator())
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
      setTimeStep();
      // compute turbulent viscosity
      _turbViscIterator.iterate();
      // TODO WS2: communicate turbulent viscosity values
      // set global boundary values for the turbulent viscosity
      _wallTurbVisIterator.iterate();
      // compute fgh for turbulent case
      _fghTurbIterator.iterate();
      // set global boundary values
      _wallFGHIterator.iterate();
      // compute the right hand side
      _rhsIterator.iterate();
      // solve for pressure
      _solver.solve();
      // TODO WS2: communicate pressure values
      // compute velocity
      _velocityIterator.iterate();
      // set obstacle boundaries
      _obstacleIterator.iterate();
      // TODO WS2: communicate velocity values
      // Iterate for velocities on the boundary
      _wallVelocityIterator.iterate();
    }

  protected:
    virtual void setTimeStep(){
      // TODO
      // iterate stencil (e.g. MinDtStencil) over all cells to find smallest dt from formula f
      // f: equation (12) from work sheet p.8, where Re=1/(nu+nuT)
      // then communicate time step to all ranks

      FLOAT localMin, globalMin;

      _minDtStencil.reset();
      _minDtIterator.iterate();

      localMin = _minDtStencil.getMinValue();

      globalMin = MY_FLOAT_MAX;
      MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

      _parameters.timestep.dt = globalMin;
      _parameters.timestep.dt *= _parameters.timestep.tau;
    }

    GlobalBoundaryIterator<TurbulentFlowField> createGlobalBoundaryTurbViscIterator(){
      BoundaryStencil<TurbulentFlowField> * stencils[6];
      // if (_parameters.simulation.scenario == "channel"){
        stencils[0] = new BFInputTurbViscosityStencil(_parameters);
        stencils[1] = new NeumannTurbViscosityBoundaryStencil(_parameters);
        stencils[2] = new MovingWallTurbViscosityStencil(_parameters);
        stencils[3] = stencils[2];
        stencils[4] = stencils[2];
        stencils[5] = stencils[2];
      // }

      if (_parameters.geometry.dim == 2){
        return GlobalBoundaryIterator<TurbulentFlowField>(_turbFlowField, _parameters,
                *(stencils[0]),*(stencils[1]),*(stencils[2]),*(stencils[3]));
      }else{
        return GlobalBoundaryIterator<TurbulentFlowField>(_turbFlowField, _parameters,
                *(stencils[0]),*(stencils[1]),*(stencils[2]),*(stencils[3]),*(stencils[4]),*(stencils[5]));
      }
    }
};

#endif // _TURBULENT_SIMULATION_H_
