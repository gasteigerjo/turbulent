#### The petsc environment can also be fixed in the makefile
# Using the externally defined $PETSC_DIR:
include ${PETSC_DIR}/conf/variables
# For MAC-Cluster:
# include ${PETSC_DIR}/conf/petscvariables
# for Thomas:
#PETSC_DIR = /opt/petsc/petsc-3.5.4
# PETSC_ARCH =

# default gnu compiler (currently not used)
# CC = g++
# compiler wrapper for mac-cluster
# CC = mpiCC
#CFLAGS = -Wall -Werror -O3 -xHost -unroll
# compiler on Ubuntu
CC = mpic++
CFLAGS = -Wall -O3 -Wno-unknown-pragmas -Werror
SRCDIR = ./
INCLUDE = -I. -Istencils ${PETSC_CC_INCLUDES}



NSMAIN = main.o

OBJ = DataStructures.o Configuration.o 3rdparty/tinyxml2/tinyxml2.o SimpleTimer.o

NSOBJ = FlowField.o LinearSolver.o Meshsize.o\
stencils/MaxUStencil.o stencils/MovingWallStencils.o stencils/PeriodicBoundaryStencils.o\
stencils/FGHStencil.o solvers/SORSolver.o solvers/PetscSolver.o \
stencils/RHSStencil.o stencils/VelocityStencil.o \
stencils/PressureBufferFillStencil.o stencils/PressureBufferReadStencil.o\
stencils/VelocityBufferFillStencil.o stencils/VelocityBufferReadStencil.o\
parallelManagers/PetscParallelConfiguration.o\
parallelManagers/PetscParallelManager.o\
GlobalBoundaryFactory.o\
stencils/BFStepInitStencil.o stencils/NeumannBoundaryStencils.o stencils/BFInputStencils.o stencils/ObstacleStencil.o\
TurbulentFlowField.o \
stencils/PostStencil.o stencils/TurbulentPostStencil.o \
VtkOutput.o \
Checkpoint.o \
stencils/FGHTurbStencil.o stencils/TurbViscosityStencil.o stencils/DistNearestWallStencil.o \
stencils/MinDtStencil.o stencils/TurbViscosityBoundaryStencil.o \
stencils/TurbViscosityBufferFillStencil.o stencils/TurbViscosityBufferReadStencil.o \
parallelManagers/PetscTurbulentParallelManager.o \

all: ns chkpt_to_vtk

ns: $(OBJ) $(NSOBJ) $(NSMAIN)
	$(CC) -o ns $(OBJ) $(NSOBJ) $(NSMAIN) $(PETSC_KSP_LIB) -lstdc++ $(CFLAGS)

chkpt_to_vtk: $(OBJ) $(NSOBJ) chkpt_to_vtk.o
	$(CC) -o chkpt_to_vtk $(OBJ) $(NSOBJ) chkpt_to_vtk.o $(PETSC_KSP_LIB) -lstdc++ $(CFLAGS) -Dchkpt_to_vtk

chkpt_to_vtk.o: main.cpp
	$(CC) -c $(CFLAGS) $(INCLUDE) -o chkpt_to_vtk.o main.cpp $(PETSC_KSP_LIB) -lstdc++ -Dchkpt_to_vtk

%.o: %.cpp
	$(CC) -c $(CFLAGS) $(INCLUDE) -o $*.o $*.cpp $(PETSC_KSP_LIB) -lstdc++

cleanall clean:
	for name in  ns main.o $(NSOBJ) $(OBJ) ; do \
	if [ -f $$name ]; then rm $$name; fi; \
	done;
