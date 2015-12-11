#include "PetscParallelManager.h"

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
            // 2D - complete edges, including the 3 ghost layers.
            // pressure points
            presBufSizeLR = Ny+3;
            presBufSizeTB = Nx+3;
            presBufSizeFB = 0;
            // velocity vectors (2D)
            velBufSizeLR  = 2*(Ny+3);
            velBufSizeTB  = 2*(Nx+3);
            velBufSizeFB  = 0;
        } else if (_parameters.geometry.dim == 3) {
            // 3D - complete planes, including the 3 ghost layers.
            // pressure points
            presBufSizeLR = (Ny+3)*(Nz+3);
            presBufSizeTB = (Nx+3)*(Nz+3);
            presBufSizeFB = (Nx+3)*(Ny+3);
            // velocity vectors (3D)
            velBufSizeLR  = 3*(Ny+3)*(Nz+3);
            velBufSizeTB  = 3*(Nx+3)*(Nz+3);
            velBufSizeFB  = 3*(Nx+3)*(Ny+3);
        }

        // Initialize the buffers
        // Pressure send buffers
        _pressuresLeftSend      = (FLOAT*) calloc(presBufSizeLR, sizeof(FLOAT));
        _pressuresRightSend     = (FLOAT*) calloc(2*presBufSizeLR, sizeof(FLOAT));
        _pressuresBottomSend    = (FLOAT*) calloc(presBufSizeTB, sizeof(FLOAT));
        _pressuresTopSend       = (FLOAT*) calloc(2*presBufSizeTB, sizeof(FLOAT));
        if (_parameters.geometry.dim == 3) {
            _pressuresFrontSend     = (FLOAT*) calloc(presBufSizeFB, sizeof(FLOAT));
            _pressuresBackSend      = (FLOAT*) calloc(2*presBufSizeFB, sizeof(FLOAT));
        }

        // Pressure receive buffers
        _pressuresLeftRecv      = (FLOAT*) calloc(2*presBufSizeLR, sizeof(FLOAT));
        _pressuresRightRecv     = (FLOAT*) calloc(presBufSizeLR, sizeof(FLOAT));
        _pressuresBottomRecv    = (FLOAT*) calloc(2*presBufSizeTB, sizeof(FLOAT));
        _pressuresTopRecv       = (FLOAT*) calloc(presBufSizeTB, sizeof(FLOAT));
        if (_parameters.geometry.dim == 3) {
            _pressuresFrontRecv     = (FLOAT*) calloc(2*presBufSizeFB, sizeof(FLOAT));
            _pressuresBackRecv      = (FLOAT*) calloc(presBufSizeFB, sizeof(FLOAT));
        }

        // Velocities send buffers
        _velocitiesLeftSend      = (FLOAT*) calloc(velBufSizeLR, sizeof(FLOAT));
        _velocitiesRightSend     = (FLOAT*) calloc(2*velBufSizeLR, sizeof(FLOAT));
        _velocitiesBottomSend    = (FLOAT*) calloc(velBufSizeTB, sizeof(FLOAT));
        _velocitiesTopSend       = (FLOAT*) calloc(2*velBufSizeTB, sizeof(FLOAT));
        if (_parameters.geometry.dim == 3) {
            _velocitiesFrontSend     = (FLOAT*) calloc(velBufSizeFB, sizeof(FLOAT));
            _velocitiesBackSend      = (FLOAT*) calloc(2*velBufSizeFB, sizeof(FLOAT));
        }

        // Velocities receive buffers
        _velocitiesLeftRecv      = (FLOAT*) calloc(2*velBufSizeLR, sizeof(FLOAT));
        _velocitiesRightRecv     = (FLOAT*) calloc(velBufSizeLR, sizeof(FLOAT));
        _velocitiesBottomRecv    = (FLOAT*) calloc(2*velBufSizeTB, sizeof(FLOAT));
        _velocitiesTopRecv       = (FLOAT*) calloc(velBufSizeTB, sizeof(FLOAT));
        if (_parameters.geometry.dim == 3) {
            _velocitiesFrontRecv     = (FLOAT*) calloc(2*velBufSizeFB, sizeof(FLOAT));
            _velocitiesBackRecv      = (FLOAT*) calloc(velBufSizeFB, sizeof(FLOAT));
        }

        int lowOffset = 0;
        int highOffset = 0;

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
        _parallelBoundaryPressureFillIterator = new ParallelBoundaryIterator<FlowField>(_flowField, _parameters, *_pressureBufferFillStencil, lowOffset, highOffset);
        _parallelBoundaryPressureReadIterator = new ParallelBoundaryIterator<FlowField>(_flowField, _parameters, *_pressureBufferReadStencil, lowOffset, highOffset);
        _parallelBoundaryVelocityFillIterator = new ParallelBoundaryIterator<FlowField>(_flowField, _parameters, *_velocityBufferFillStencil, lowOffset, highOffset);
        _parallelBoundaryVelocityReadIterator = new ParallelBoundaryIterator<FlowField>(_flowField, _parameters, *_velocityBufferReadStencil, lowOffset, highOffset);

    }


void PetscParallelManager::communicatePressure(){

    MPI_Status comm_status;
    MPI_Request send_requestL, send_requestR, send_requestBt, send_requestT, send_requestF, send_requestBk;
    MPI_Request recv_requestL, recv_requestR, recv_requestBt, recv_requestT, recv_requestF, recv_requestBk;

    // ------------------------------------------------------------------------
    // Communicate pressure in the X-dimension:
    // Fill the pressure send buffers L/R
    _parallelBoundaryPressureFillIterator->iterate(0);

    // Left --> + --> Right
    MPI_Isend(_pressuresRightSend, 2*presBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 1, PETSC_COMM_WORLD, &send_requestR);
    // printf("[rank %d] L-->R Send pressure request posted.\n", _parameters.parallel.rank);
    MPI_Irecv(_pressuresLeftRecv, 2*presBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,   1, PETSC_COMM_WORLD, &recv_requestL);
    // printf("[rank %d] L-->R Recv pressure request posted.\n", _parameters.parallel.rank);

    // Left <-- + <-- Right
    MPI_Isend(_pressuresLeftSend,  presBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  2, PETSC_COMM_WORLD, &send_requestL);
    // printf("[rank %d] L<--R Send pressure request posted.\n", _parameters.parallel.rank);
    MPI_Irecv(_pressuresRightRecv, presBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 2, PETSC_COMM_WORLD, &recv_requestR);
    // printf("[rank %d] L<--R Recv pressure request posted.\n", _parameters.parallel.rank);

    // Wait for all the requests of the current rank to finish before reading
    MPI_Wait(&recv_requestL, &comm_status);
    // printf("[rank %d] L-->R RECV pressure request COMPLETED.\n", _parameters.parallel.rank);
    MPI_Wait(&recv_requestR, &comm_status);
    // printf("[rank %d] L<--R RECV pressure request COMPLETED.\n", _parameters.parallel.rank);

    // Read the pressure receive buffers L/R
    _parallelBoundaryPressureReadIterator->iterate(0);
    // ------------------------------------------------------------------------



    // ------------------------------------------------------------------------
    // Communicate pressure in the Y-dimension:
    // Fill the pressure send buffers T/B
    _parallelBoundaryPressureFillIterator->iterate(1);

    // Bottom --> + --> Top
    MPI_Isend(_pressuresTopSend,    2*presBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    3, PETSC_COMM_WORLD, &send_requestT);
    // printf("[rank %d] B-->T Send pressure request posted.\n", _parameters.parallel.rank);
    MPI_Irecv(_pressuresBottomRecv, 2*presBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 3, PETSC_COMM_WORLD, &recv_requestBt);
    // printf("[rank %d] B-->T Recv pressure request posted.\n", _parameters.parallel.rank);

    // Bottom <-- + <-- Top
    MPI_Isend(_pressuresBottomSend, presBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 4, PETSC_COMM_WORLD, &send_requestBt);
    // printf("[rank %d] B<--T Send pressure request posted.\n", _parameters.parallel.rank);
    MPI_Irecv(_pressuresTopRecv,    presBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    4, PETSC_COMM_WORLD, &recv_requestT);
    // printf("[rank %d] B<--T Recv pressure request posted.\n", _parameters.parallel.rank);

    // Wait for all the requests of the current rank to finish before reading
    MPI_Wait(&recv_requestBt, &comm_status);
    // printf("[rank %d] B-->T RECV pressure request COMPLETED.\n", _parameters.parallel.rank);
    MPI_Wait(&recv_requestT,  &comm_status);
    // printf("[rank %d] B<--T RECV pressure request COMPLETED.\n", _parameters.parallel.rank);

    // Read the pressure receive buffers T/B
    _parallelBoundaryPressureReadIterator->iterate(1);
    // ------------------------------------------------------------------------



    // ------------------------------------------------------------------------
    // Communicate pressure in the Z-dimension  (z increases from front to back):
    if (_parameters.geometry.dim == 3) {

        // Fill the pressure send buffers F/B
        _parallelBoundaryPressureFillIterator->iterate(2);

        // Front --> + --> Back
        MPI_Isend(_pressuresBackSend,  2*presBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  5, PETSC_COMM_WORLD, &send_requestBk);
        // printf("[rank %d] F-->B Send pressure request posted.\n", _parameters.parallel.rank);
        MPI_Irecv(_pressuresFrontRecv, 2*presBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 5, PETSC_COMM_WORLD, &recv_requestF);
        // printf("[rank %d] F-->B Recv pressure request posted.\n", _parameters.parallel.rank);

        // Front <-- + <-- Back
        MPI_Isend(_pressuresFrontSend, presBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 6, PETSC_COMM_WORLD, &send_requestF);
        // printf("[rank %d] F<--B Send pressure request posted.\n", _parameters.parallel.rank);
        MPI_Irecv(_pressuresBackRecv,  presBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  6, PETSC_COMM_WORLD, &recv_requestBk);
        // printf("[rank %d] F<--B Recv pressure request posted.\n", _parameters.parallel.rank);

        // Wait for all the requests of the current rank to finish before reading
        MPI_Wait(&recv_requestF,  &comm_status);
        // printf("[rank %d] F-->B RECV pressure request COMPLETED.\n", _parameters.parallel.rank);
        MPI_Wait(&recv_requestBk, &comm_status);
        // printf("[rank %d] F<--B RECV pressure request COMPLETED.\n", _parameters.parallel.rank);

        // Read the pressure receive buffers F/B
        _parallelBoundaryPressureReadIterator->iterate(2);
    }
    // ------------------------------------------------------------------------

    // Wait for all the sends of this rank to complete
    MPI_Wait(&send_requestL,  &comm_status);
    MPI_Wait(&send_requestR,  &comm_status);
    MPI_Wait(&send_requestBt, &comm_status);
    MPI_Wait(&send_requestT,  &comm_status);
    MPI_Wait(&send_requestF,  &comm_status);
    MPI_Wait(&send_requestBk, &comm_status);

}

void PetscParallelManager::communicateVelocities(){

    MPI_Status comm_status;
    MPI_Request send_requestL, send_requestR, send_requestBt, send_requestT, send_requestF, send_requestBk;
    MPI_Request recv_requestL, recv_requestR, recv_requestBt, recv_requestT, recv_requestF, recv_requestBk;

    // ------------------------------------------------------------------------
    // Communicate velocities in the X-dimension:
    // Fill the velocity send buffers L/R
    _parallelBoundaryVelocityFillIterator->iterate(0);

    // Left --> + --> Right
    MPI_Isend(_velocitiesRightSend, 2*velBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 11, PETSC_COMM_WORLD, &send_requestR);
    MPI_Irecv(_velocitiesLeftRecv,  2*velBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  11, PETSC_COMM_WORLD, &recv_requestL);

    // Left <-- + <-- Right
    MPI_Isend(_velocitiesLeftSend,  velBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  12, PETSC_COMM_WORLD, &send_requestL);
    MPI_Irecv(_velocitiesRightRecv, velBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 12, PETSC_COMM_WORLD, &recv_requestR);

    // Wait for all the requests of the current rank to finish before reading
    MPI_Wait(&recv_requestL, &comm_status);
    MPI_Wait(&recv_requestR, &comm_status);

    // Read the velocity receive buffers L/R
    _parallelBoundaryVelocityReadIterator->iterate(0);
    // ------------------------------------------------------------------------



    // ------------------------------------------------------------------------
    // Communicate velocities in the Y-dimension:
    // Fill the velocity send buffers T/B
    _parallelBoundaryVelocityFillIterator->iterate(1);

    // Bottom --> + --> Top
    MPI_Isend(_velocitiesTopSend,    2*velBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    13, PETSC_COMM_WORLD, &send_requestT);
    MPI_Irecv(_velocitiesBottomRecv, 2*velBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 13, PETSC_COMM_WORLD, &recv_requestBt);

    // Bottom <-- + <-- Top
    MPI_Isend(_velocitiesBottomSend, velBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 14, PETSC_COMM_WORLD, &send_requestBt);
    MPI_Irecv(_velocitiesTopRecv,    velBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    14, PETSC_COMM_WORLD, &recv_requestT);

    // Wait for all the requests of the current rank to finish before reading
    MPI_Wait(&recv_requestBt, &comm_status);
    MPI_Wait(&recv_requestT,  &comm_status);

    // Read the velocity receive buffers T/B
    _parallelBoundaryVelocityReadIterator->iterate(1);
    // ------------------------------------------------------------------------



    // ------------------------------------------------------------------------
    // Communicate velocities in the Z-dimension (z increases from front to back):
    if (_parameters.geometry.dim == 3) {

        // Fill the pressure send buffers F/B
        _parallelBoundaryVelocityFillIterator->iterate(2);

        // Front --> + --> Back
        MPI_Isend(_velocitiesBackSend,  2*velBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  15, PETSC_COMM_WORLD, &send_requestBk);
        MPI_Irecv(_velocitiesFrontRecv, 2*velBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 15, PETSC_COMM_WORLD, &recv_requestF);

        // Front <-- + <-- Back
        MPI_Isend(_velocitiesFrontSend, velBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 16, PETSC_COMM_WORLD, &send_requestF);
        MPI_Irecv(_velocitiesBackRecv,  velBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  16, PETSC_COMM_WORLD, &recv_requestBk);

        // Wait for all the requests of the current rank to finish before reading
        MPI_Wait(&recv_requestF,  &comm_status);
        MPI_Wait(&recv_requestBk, &comm_status);

        // Read the velocity receive buffers F/B
        _parallelBoundaryVelocityReadIterator->iterate(2);
    }
    // ------------------------------------------------------------------------

    // Wait for all the sends of this rank to complete
    MPI_Wait(&send_requestL,  &comm_status);
    MPI_Wait(&send_requestR,  &comm_status);
    MPI_Wait(&send_requestBt, &comm_status);
    MPI_Wait(&send_requestT,  &comm_status);
    MPI_Wait(&send_requestF,  &comm_status);
    MPI_Wait(&send_requestBk, &comm_status);

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
