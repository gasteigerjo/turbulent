#include "PetscParallelManager.h"

// TODO: Open issues and notes
// - Are the offsets set correctly? We really need to check this!
// - We needed to initialize the stencils and iterators before the body of the constructor. Do some discussion.
// - Is there any problem with the calloc-free method, instead of the new-delete?

// Constructor
PetscParallelManager::PetscParallelManager(const Parameters & parameters, FlowField & flowField):
_parameters(parameters),
_flowField(flowField)
{
        // Get the domain size
        int Nx = _parameters.parallel.localSize[0];
        int Ny = _parameters.parallel.localSize[1];
        int Nz = _parameters.parallel.localSize[2];

        // Define the buffer sizes
        if (_parameters.geometry.dim == 2) {
            // 2D - complete edges, including ghost layers.
            presBufSizeLR = Ny+2;
            presBufSizeTB = Nx+2;
            presBufSizeFB = 0;
            velBufSizeLR  = 2*(Ny+2);
            velBufSizeTB  = 2*(Nx+2);
            velBufSizeFB  = 0;
        } else if (_parameters.geometry.dim == 3) {
            // 3D - complete planes, including ghost layers.
            presBufSizeLR = (Ny+2)*(Nz+2);
            presBufSizeTB = (Nx+2)*(Nz+2);
            presBufSizeFB = (Nx+2)*(Ny+2);
            velBufSizeLR  = 3*(Ny+2)*(Nz+2);
            velBufSizeTB  = 3*(Nx+2)*(Nz+2);
            velBufSizeFB  = 3*(Nx+2)*(Ny+2);
        }

        // Initialize the buffers
        // Pressure send buffers
        _pressuresLeftSend      = (FLOAT*) calloc(presBufSizeLR, sizeof(FLOAT));
        _pressuresRightSend     = (FLOAT*) calloc(presBufSizeLR, sizeof(FLOAT));
        _pressuresBottomSend    = (FLOAT*) calloc(presBufSizeTB, sizeof(FLOAT));
        _pressuresTopSend       = (FLOAT*) calloc(presBufSizeTB, sizeof(FLOAT));
        if (_parameters.geometry.dim == 3) {
            _pressuresFrontSend     = (FLOAT*) calloc(presBufSizeFB, sizeof(FLOAT));
            _pressuresBackSend      = (FLOAT*) calloc(presBufSizeFB, sizeof(FLOAT));
        }

        // Pressure receive buffers
        _pressuresLeftRecv      = (FLOAT*) calloc(presBufSizeLR, sizeof(FLOAT));
        _pressuresRightRecv     = (FLOAT*) calloc(presBufSizeLR, sizeof(FLOAT));
        _pressuresBottomRecv    = (FLOAT*) calloc(presBufSizeTB, sizeof(FLOAT));
        _pressuresTopRecv       = (FLOAT*) calloc(presBufSizeTB, sizeof(FLOAT));
        if (_parameters.geometry.dim == 3) {
            _pressuresFrontRecv     = (FLOAT*) calloc(presBufSizeFB, sizeof(FLOAT));
            _pressuresBackRecv      = (FLOAT*) calloc(presBufSizeFB, sizeof(FLOAT));
        }

        // Velocities send buffers
        _velocitiesLeftSend      = (FLOAT*) calloc(velBufSizeLR, sizeof(FLOAT));
        _velocitiesRightSend     = (FLOAT*) calloc(velBufSizeLR, sizeof(FLOAT));
        _velocitiesBottomSend    = (FLOAT*) calloc(velBufSizeTB, sizeof(FLOAT));
        _velocitiesTopSend       = (FLOAT*) calloc(velBufSizeTB, sizeof(FLOAT));
        if (_parameters.geometry.dim == 3) {
            _velocitiesFrontSend     = (FLOAT*) calloc(velBufSizeFB, sizeof(FLOAT));
            _velocitiesBackSend      = (FLOAT*) calloc(velBufSizeFB, sizeof(FLOAT));
        }

        // Velocities receive buffers
        _velocitiesLeftRecv      = (FLOAT*) calloc(velBufSizeLR, sizeof(FLOAT));
        _velocitiesRightRecv     = (FLOAT*) calloc(velBufSizeLR, sizeof(FLOAT));
        _velocitiesBottomRecv    = (FLOAT*) calloc(velBufSizeTB, sizeof(FLOAT));
        _velocitiesTopRecv       = (FLOAT*) calloc(velBufSizeTB, sizeof(FLOAT));
        if (_parameters.geometry.dim == 3) {
            _velocitiesFrontRecv     = (FLOAT*) calloc(velBufSizeFB, sizeof(FLOAT));
            _velocitiesBackRecv      = (FLOAT*) calloc(velBufSizeFB, sizeof(FLOAT));
        }

        int lowOffset = 2;
        int highOffset = -1;
        // Construct the stencils
        if (_parameters.geometry.dim == 2) {
            _pressureBufferFillStencil = new PressureBufferFillStencil(_parameters, _pressuresLeftSend,  _pressuresRightSend,  _pressuresBottomSend,  _pressuresTopSend, lowOffset);
            _pressureBufferReadStencil = new PressureBufferReadStencil(_parameters, _pressuresLeftRecv,  _pressuresRightRecv,  _pressuresBottomRecv,  _pressuresTopRecv, lowOffset);
            _velocityBufferFillStencil = new VelocityBufferFillStencil(_parameters, _velocitiesLeftSend, _velocitiesRightSend, _velocitiesBottomSend, _velocitiesTopSend, lowOffset);
            _velocityBufferReadStencil = new VelocityBufferReadStencil(_parameters, _velocitiesLeftRecv, _velocitiesRightRecv, _velocitiesBottomRecv, _velocitiesTopRecv, lowOffset);
        } else if (_parameters.geometry.dim == 3) {
            _pressureBufferFillStencil = new PressureBufferFillStencil(_parameters, _pressuresLeftSend,  _pressuresRightSend,  _pressuresBottomSend,  _pressuresTopSend,  _pressuresFrontSend,  _pressuresBackSend, lowOffset);
            _pressureBufferReadStencil = new PressureBufferReadStencil(_parameters, _pressuresLeftRecv,  _pressuresRightRecv,  _pressuresBottomRecv,  _pressuresTopRecv,  _pressuresFrontRecv,  _pressuresBackRecv, lowOffset);
            _velocityBufferFillStencil = new VelocityBufferFillStencil(_parameters, _velocitiesLeftSend, _velocitiesRightSend, _velocitiesBottomSend, _velocitiesTopSend, _velocitiesFrontSend, _velocitiesBackSend, lowOffset);
            _velocityBufferReadStencil = new VelocityBufferReadStencil(_parameters, _velocitiesLeftRecv, _velocitiesRightRecv, _velocitiesBottomRecv, _velocitiesTopRecv, _velocitiesFrontRecv, _velocitiesBackRecv, lowOffset);
        }

        // Construct the iterators
        _parallelBoundaryPressureFillIterator = new ParallelBoundaryIterator(_flowField, _parameters, _pressureBufferFillStencil, lowOffset, highOffset);
        _parallelBoundaryPressureReadIterator = new ParallelBoundaryIterator(_flowField, _parameters, _pressureBufferReadStencil, lowOffset, highOffset);
        _parallelBoundaryVelocityFillIterator = new ParallelBoundaryIterator(_flowField, _parameters, _velocityBufferFillStencil, lowOffset, highOffset);
        _parallelBoundaryVelocityReadIterator = new ParallelBoundaryIterator(_flowField, _parameters, _velocityBufferReadStencil, lowOffset, highOffset);

    }

void PetscParallelManager::communicatePressure(){

    // Fill the pressure send buffers
    printf("BoundaryPressureFillIterator... (rank %d)\n", _parameters.parallel.rank);
    _parallelBoundaryPressureFillIterator->iterate();

    MPI_Status comm_status;

    // Communicate in the X-dimension:
    // Left --> + --> Right
    printf("Communicating L-->R... (rank %d)\n", _parameters.parallel.rank);
    MPI_Sendrecv(_pressuresRightSend, presBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 1,
                 _pressuresLeftRecv,  presBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  1,
                 PETSC_COMM_WORLD, &comm_status);
    // Left <-- + <-- Right
    printf("Communicating L<--R... (rank %d)\n", _parameters.parallel.rank);
    MPI_Sendrecv(_pressuresLeftSend,  presBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  1,
                 _pressuresRightRecv, presBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 1,
                 PETSC_COMM_WORLD, &comm_status);

    // Communicate in the Y-dimension:
    // Top --> + --> Bottom
    printf("Communicating T-->B... (rank %d)\n", _parameters.parallel.rank);
    MPI_Sendrecv(_pressuresBottomSend, presBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 1,
                 _pressuresTopRecv,    presBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    1,
                 PETSC_COMM_WORLD, &comm_status);
    // Top <-- + <-- Bottom
    printf("Communicating T<--B... (rank %d)\n", _parameters.parallel.rank);
    MPI_Sendrecv(_pressuresTopSend,    presBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    1,
                 _pressuresBottomRecv, presBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 1,
                 PETSC_COMM_WORLD, &comm_status);

    if (_parameters.geometry.dim == 3) {
        // Communicate in the Z-dimension:
        // Front --> + --> Back
        printf("Communicating F-->B... (rank %d)\n", _parameters.parallel.rank);
        MPI_Sendrecv(_pressuresBackSend,  presBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  1,
                     _pressuresFrontRecv, presBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 1,
                     PETSC_COMM_WORLD, &comm_status);
        // Front <-- + <-- Back
        printf("Communicating F<--B... (rank %d)\n", _parameters.parallel.rank);
        MPI_Sendrecv(_pressuresFrontSend, presBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 1,
                     _pressuresBackRecv,  presBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  1,
                     PETSC_COMM_WORLD, &comm_status);
    }

    // Read the pressure receive buffers
    printf("BoundaryPressureReadIterator... (rank %d)\n", _parameters.parallel.rank);
    _parallelBoundaryPressureReadIterator->iterate();

}

void PetscParallelManager::communicateVelocities(){

        // Fill the velocity send buffers
        _parallelBoundaryVelocityFillIterator->iterate();

        MPI_Status comm_status;

        // Communicate in the X-dimension:
        // Left --> + --> Right
        MPI_Sendrecv(_velocitiesRightSend, velBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 1,
                     _velocitiesLeftRecv,  velBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  1,
                     PETSC_COMM_WORLD, &comm_status);
        // Left <-- + <-- Right
        MPI_Sendrecv(_velocitiesLeftSend,  velBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  1,
                     _velocitiesRightRecv, velBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 1,
                     PETSC_COMM_WORLD, &comm_status);

        // Communicate in the Y-dimension:
        // Top --> + --> Bottom
        MPI_Sendrecv(_velocitiesBottomSend, velBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 1,
                     _velocitiesTopRecv,    velBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    1,
                     PETSC_COMM_WORLD, &comm_status);
        // Top <-- + <-- Bottom
        MPI_Sendrecv(_velocitiesTopSend,    velBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    1,
                     _velocitiesBottomRecv, velBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 1,
                     PETSC_COMM_WORLD, &comm_status);

        if (_parameters.geometry.dim == 3) {
            // Communicate in the Z-dimension:
            // Front --> + --> Back
            MPI_Sendrecv(_velocitiesBackSend,  velBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  1,
                         _velocitiesFrontRecv, velBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 1,
                         PETSC_COMM_WORLD, &comm_status);
            // Front <-- + <-- Back
            MPI_Sendrecv(_velocitiesFrontSend, velBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 1,
                         _velocitiesBackRecv,  velBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  1,
                         PETSC_COMM_WORLD, &comm_status);
        }

        // Read the velocity receive buffers
        _parallelBoundaryVelocityReadIterator->iterate();

}

// Destructor
PetscParallelManager::~PetscParallelManager() {
    // Free the pressure buffers
    free(_pressuresLeftSend);
    free(_pressuresRightSend);
    free(_pressuresBottomSend);
    free(_pressuresTopSend);

    free(_pressuresLeftRecv);
    free(_pressuresRightRecv);
    free(_pressuresBottomRecv);
    free(_pressuresTopRecv);

    if (_parameters.geometry.dim == 3) {
        free(_pressuresFrontSend);
        free(_pressuresBackSend);
        free(_pressuresFrontRecv);
        free(_pressuresBackRecv);
    }

    // Free the velocities buffers
    free(_velocitiesLeftSend);
    free(_velocitiesRightSend);
    free(_velocitiesBottomSend);
    free(_velocitiesTopSend);

    free(_velocitiesLeftRecv);
    free(_velocitiesRightRecv);
    free(_velocitiesBottomRecv);
    free(_velocitiesTopRecv);

    if (_parameters.geometry.dim == 3) {
        free(_velocitiesFrontSend);
        free(_velocitiesBackSend);
        free(_velocitiesFrontRecv);
        free(_velocitiesBackRecv);
    }
}
