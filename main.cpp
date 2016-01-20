#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "Configuration.h"
#include "Simulation.h"
#include "parallelManagers/PetscParallelConfiguration.h"
#include "MeshsizeFactory.h"
#include <iomanip>
#include "TurbulentSimulation.h"
#include "SimpleTimer.h"

int main (int argc, char *argv[]) {

    // Parallelization related. Initialize and identify
    // ---------------------------------------------------
    int rank;   // This processor's identifier
    int nproc;  // Number of processors in the group
    PetscInitialize(&argc, &argv, "petsc_commandline_arg", PETSC_NULL);
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    std::cout << "Rank: " << rank << ", Nproc: " << nproc << std::endl;
    //----------------------------------------------------


    // read configuration and store information in parameters object
    Configuration configuration(argv[1]);
    Parameters parameters;
    configuration.loadParameters(parameters);
    PetscParallelConfiguration parallelConfiguration(parameters);
    MeshsizeFactory::getInstance().initMeshsize(parameters);
    FlowField *flowField = NULL;
    Simulation *simulation = NULL;
    SimpleTimer timer = SimpleTimer();

    #ifdef DEBUG
    std::cout << "Processor " << parameters.parallel.rank << " with index ";
    std::cout << parameters.parallel.indices[0] << ",";
    std::cout << parameters.parallel.indices[1] << ",";
    std::cout << parameters.parallel.indices[2];
    std::cout <<    " is computing the size of its subdomain and obtains ";
    std::cout << parameters.parallel.localSize[0] << ", ";
    std::cout << parameters.parallel.localSize[1] << " and ";
    std::cout << parameters.parallel.localSize[2];
    std::cout << ". Left neighbour: " << parameters.parallel.leftNb;
    std::cout << ", right neighbour: " << parameters.parallel.rightNb;
    std::cout << std::endl;
    std::cout << "Min. meshsizes: " << parameters.meshsize->getDxMin() << ", " << parameters.meshsize->getDyMin() << ", " << parameters.meshsize->getDzMin() << std::endl;
    #endif

    // DEBUG
    std::cout << "Checkpoint iterations: " << parameters.checkpoint.iterations << ", directory: " << parameters.checkpoint.directory << ", prefix: " << parameters.checkpoint.prefix << ", cleanDirectory:" << parameters.checkpoint.cleanDirectory << std::endl;
    std::cout << "Restart filename: " << parameters.restart.filename << std::endl;

    // initialise simulation
    if (parameters.simulation.type=="turbulence"){
      if(rank==0){ std::cout << "Start RANS simulation in " << parameters.geometry.dim << "D" << std::endl; }
      // WS2: initialise turbulent flow field and turbulent simulation object
      TurbulentFlowField *turbFlowField = NULL;
      turbFlowField = new TurbulentFlowField(parameters);
      flowField = turbFlowField;
      if(flowField == NULL){ handleError(1, "flowField==NULL!"); }
      simulation = new TurbulentSimulation(parameters,*turbFlowField);
      // handleError(1,"Turbulence currently not supported yet!");
    } else if (parameters.simulation.type=="dns"){
      if(rank==0){ std::cout << "Start DNS simulation in " << parameters.geometry.dim << "D" << std::endl; }
      flowField = new FlowField(parameters);
      if(flowField == NULL){ handleError(1, "flowField==NULL!"); }
      simulation = new Simulation(parameters,*flowField);
    } else {
      handleError(1, "Unknown simulation type! Currently supported: dns, turbulence");
    }
    // call initialization of simulation (initialize flow field)
    if(simulation == NULL){ handleError(1, "simulation==NULL!"); }
    simulation->initializeFlowField();
    //flowField->getFlags().show();

    int timeSteps = 0;
    FLOAT time = 0.0;

    // Read the restart data
    if(parameters.restart.filename != "") {
        simulation->readCheckpoint(timeSteps, time);
        printf(" ++++ timestep: %d, time: %f\n", timeSteps, time); //DEBUG
    }

    FLOAT lastPlotTime = time;
    int lastCheckpointIter = timeSteps;
    FLOAT timeStdOut=parameters.stdOut.interval;

    // WS1: plot initial state
    simulation->plotVTK(0);

    // initialize the region timers
    FLOAT time_loop  = 0; FLOAT time_loop_tot  = 0;
    FLOAT time_solve = 0; FLOAT time_solve_tot = 0;
    FLOAT time_comm  = 0; FLOAT time_comm_tot  = 0;

    // clean the checkpoints directory if needed
    if (parameters.checkpoint.cleanDirectory) {
        simulation->cleandirCheckpoint();
    }

    // create the first checkpoint
    // simulation->createCheckpoint(timeSteps, time);

    // start the global timer
    timer.start();

    // time loop
    while (time < parameters.simulation.finalTime){

      simulation->solveTimestep(time_solve, time_comm);

      time += parameters.timestep.dt;
      timeSteps++;

      // std-out: terminal info
      if ( (rank==0) && (timeStdOut <= time) ){
          std::cout << "Current time: " << time << "\ttimestep: " <<
                        parameters.timestep.dt << "\titeration: " << timeSteps <<std::endl << std::endl;
          timeStdOut += parameters.stdOut.interval;
      }

      // DEBUG: Currently, restarting the simulation will overwrite the last checkpoint. Change that!
      if (lastCheckpointIter + parameters.checkpoint.iterations <= timeSteps) {
          simulation->createCheckpoint(timeSteps, time);
          lastCheckpointIter += parameters.checkpoint.iterations;
      }

      // WS1: trigger VTK output
      if (lastPlotTime + parameters.vtk.interval <= time) {
          simulation->plotVTK(timeSteps); // TODO Change to time?
          lastPlotTime += parameters.vtk.interval;
      }
    }

    // take computation time
    time_loop = timer.getTimeAndContinue();
    printf("[Rank %d] Timers (s):\tloop: %f | solve: %f  comm: %f  other: %f\n", rank, time_loop, time_solve, time_comm, time_loop-time_solve-time_comm);
    MPI_Reduce(&time_loop,  &time_loop_tot,  1, MY_MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);
    MPI_Reduce(&time_solve, &time_solve_tot, 1, MY_MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);
    MPI_Reduce(&time_comm,  &time_comm_tot,  1, MY_MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank==0) {
        printf("[Average] Timers (s):\tloop: %f | solve: %f  comm: %f  other: %f\n", time_loop_tot/nproc, time_solve_tot/nproc, time_comm_tot/nproc, (time_loop_tot-time_solve_tot-time_comm_tot)/nproc);
        std::cerr << parameters.parallel.numProcessors[0] << "x" << parameters.parallel.numProcessors[1] << "x" << parameters.parallel.numProcessors[2] << ": " << time_loop_tot/nproc << std::endl; // Output time in cerr for easy redirection into file
    }

    // WS1: plot final output
    simulation->plotVTK(timeSteps);

    delete simulation; simulation=NULL;
    delete flowField;  flowField= NULL;

    PetscFinalize();
}
