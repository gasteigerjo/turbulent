#### The petsc environment can also be fixed in the makefile
# PETSC_DIR =
# for Thomas:
PETSC_DIR = /opt/petsc/petsc-3.5.4
# PETSC_ARCH =
include ${PETSC_DIR}/conf/variables

# default gnu compiler (currently not used)
# CC = g++
# compiler wrapper for mac-cluster
# CC = mpiCC
# compiler on Ubuntu
#
#
CC = mpic++
CFLAGS = -Wall -Werror -O3 -Wno-unknown-pragmas -Wno-unused-value
SRCDIR = ./
INCLUDE = -I. -Istencils ${PETSC_CC_INCLUDES}


NSMAIN = main.o

OBJ = DataStructures.o Configuration.o 3rdparty/tinyxml2/tinyxml2.o SimpleTimer.o

NSOBJ = FlowField.o LinearSolver.o Meshsize.o\
stencils/MaxUStencil.o stencils/MovingWallStencils.o stencils/PeriodicBoundaryStencils.o\
stencils/FGHStencil.o solvers/SORSolver.o solvers/PetscSolver.o \
stencils/RHSStencil.o stencils/VelocityStencil.o \
stencils/VTKStencil.o \
stencils/PressureBufferFillStencil.o stencils/PressureBufferReadStencil.o\
stencils/VelocityBufferFillStencil.o stencils/VelocityBufferReadStencil.o\
parallelManagers/PetscParallelConfiguration.o\
parallelManagers/PetscParallelManager.o\
GlobalBoundaryFactory.o\
stencils/BFStepInitStencil.o stencils/NeumannBoundaryStencils.o stencils/BFInputStencils.o stencils/ObstacleStencil.o\
TurbulentFlowField.o \
stencils/FGHTurbStencil.o stencils/TurbViscosityStencil.o stencils/DistNearestWallStencil.o \
stencils/MinDtStencil.o stencils/TurbViscosityBoundaryStencil.o \

all: ns

ns: $(OBJ) $(NSOBJ) $(NSMAIN)
	$(CC) -o ns $(OBJ) $(NSOBJ) $(NSMAIN) $(PETSC_KSP_LIB) -lstdc++ $(CFLAGS)


%.o: %.cpp
	$(CC) -c $(CFLAGS) $(INCLUDE) -o $*.o $*.cpp $(PETSC_KSP_LIB) -lstdc++

cleanall clean:
	for name in  ns main.o $(NSOBJ) $(OBJ) ; do \
	if [ -f $$name ]; then rm $$name; fi; \
	done;
