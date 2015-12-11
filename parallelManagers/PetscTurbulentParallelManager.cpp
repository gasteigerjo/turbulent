#include "PetscTurbulentParallelManager.h"

// Constructor
PetscTurbulentParallelManager::PetscTurbulentParallelManager(const Parameters & parameters, TurbulentFlowField & turbFlowField):
// PetscParallelManager(parameters,turbFlowField),
_parameters(parameters),
_turbFlowField(turbFlowField)
{
        // Get the domain size
        int Nx = _parameters.parallel.localSize[0];
        int Ny = _parameters.parallel.localSize[1];
        int Nz = _parameters.parallel.localSize[2];

        // Define the buffer sizes
        if (_parameters.geometry.dim == 2) {
            // 2D - complete edges, including the 3 ghost layers.
            // turbulent viscosity points
            turbViscBufSizeLR = 2*(Ny+3);
            turbViscBufSizeTB = 2*(Nx+3);
            turbViscBufSizeFB = 0;
        } else if (_parameters.geometry.dim == 3) {
            // 3D - complete planes, including the 3 ghost layers.
            // turbulent viscosity points
            turbViscBufSizeLR = 2*(Ny+3)*(Nz+3);
            turbViscBufSizeTB = 2*(Nx+3)*(Nz+3);
            turbViscBufSizeFB = 2*(Nx+3)*(Ny+3);
        }

        // Initialize the buffers
        // Turbulent viscosity send buffers
        _turbViscLeftSend      = (FLOAT*) calloc(turbViscBufSizeLR, sizeof(FLOAT));
        _turbViscRightSend     = (FLOAT*) calloc(turbViscBufSizeLR, sizeof(FLOAT));
        _turbViscBottomSend    = (FLOAT*) calloc(turbViscBufSizeTB, sizeof(FLOAT));
        _turbViscTopSend       = (FLOAT*) calloc(turbViscBufSizeTB, sizeof(FLOAT));
        if (_parameters.geometry.dim == 3) {
            _turbViscFrontSend     = (FLOAT*) calloc(turbViscBufSizeFB, sizeof(FLOAT));
            _turbViscBackSend      = (FLOAT*) calloc(turbViscBufSizeFB, sizeof(FLOAT));
        }

        // Turbulent viscosity receive buffers
        _turbViscLeftRecv      = (FLOAT*) calloc(turbViscBufSizeLR, sizeof(FLOAT));
        _turbViscRightRecv     = (FLOAT*) calloc(turbViscBufSizeLR, sizeof(FLOAT));
        _turbViscBottomRecv    = (FLOAT*) calloc(turbViscBufSizeTB, sizeof(FLOAT));
        _turbViscTopRecv       = (FLOAT*) calloc(turbViscBufSizeTB, sizeof(FLOAT));
        if (_parameters.geometry.dim == 3) {
            _turbViscFrontRecv     = (FLOAT*) calloc(turbViscBufSizeFB, sizeof(FLOAT));
            _turbViscBackRecv      = (FLOAT*) calloc(turbViscBufSizeFB, sizeof(FLOAT));
        }


        int lowOffset = 0;
        int highOffset = 0;

        // Construct the stencils
        if (_parameters.geometry.dim == 2) {
            _turbViscBufferFillStencil = new TurbViscosityBufferFillStencil(_parameters, _turbViscLeftSend,  _turbViscRightSend,  _turbViscBottomSend,  _turbViscTopSend, lowOffset);
            _turbViscBufferReadStencil = new TurbViscosityBufferReadStencil(_parameters, _turbViscLeftRecv,  _turbViscRightRecv,  _turbViscBottomRecv,  _turbViscTopRecv, lowOffset);
        } else if (_parameters.geometry.dim == 3) {
            _turbViscBufferFillStencil = new TurbViscosityBufferFillStencil(_parameters, _turbViscLeftSend,  _turbViscRightSend,  _turbViscBottomSend,  _turbViscTopSend,  _turbViscFrontSend,  _turbViscBackSend, lowOffset);
            _turbViscBufferReadStencil = new TurbViscosityBufferReadStencil(_parameters, _turbViscLeftRecv,  _turbViscRightRecv,  _turbViscBottomRecv,  _turbViscTopRecv,  _turbViscFrontRecv,  _turbViscBackRecv, lowOffset);
        }

        // Construct the iterators
        _parallelBoundaryTurbViscFillIterator = new ParallelBoundaryIterator<TurbulentFlowField>(_turbFlowField, _parameters, *_turbViscBufferFillStencil, lowOffset, highOffset);
        _parallelBoundaryTurbViscReadIterator = new ParallelBoundaryIterator<TurbulentFlowField>(_turbFlowField, _parameters, *_turbViscBufferReadStencil, lowOffset, highOffset);

    }


void PetscTurbulentParallelManager::communicateTurbViscosity(){

    MPI_Status comm_status;
    MPI_Request send_requestL, send_requestR, send_requestBt, send_requestT, send_requestF, send_requestBk;
    MPI_Request recv_requestL, recv_requestR, recv_requestBt, recv_requestT, recv_requestF, recv_requestBk;

    // ------------------------------------------------------------------------
    // Communicate in the X-dimension:
    // Fill the turbulent viscosity send buffers L/R
    _parallelBoundaryTurbViscFillIterator->iterate(0);

    // Left --> + --> Right
    MPI_Isend(_turbViscRightSend, turbViscBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 1, PETSC_COMM_WORLD, &send_requestR);
    MPI_Irecv(_turbViscLeftRecv,  turbViscBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  1, PETSC_COMM_WORLD, &recv_requestL);

    // Left <-- + <-- Right
    MPI_Isend(_turbViscLeftSend,  turbViscBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  2, PETSC_COMM_WORLD, &send_requestL);
    MPI_Irecv(_turbViscRightRecv, turbViscBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 2, PETSC_COMM_WORLD, &recv_requestR);

    // Wait for all the requests of the current rank to finish before reading
    MPI_Wait(&recv_requestL, &comm_status);
    MPI_Wait(&recv_requestR, &comm_status);

    // Read the turbulent viscosity receive buffers L/R
    _parallelBoundaryTurbViscReadIterator->iterate(0);
    // ------------------------------------------------------------------------



    // ------------------------------------------------------------------------
    // Communicate in the Y-dimension:
    // Fill the turbulent viscosity send buffers T/B
    _parallelBoundaryTurbViscFillIterator->iterate(1);

    // Bottom --> + --> Top
    MPI_Isend(_turbViscTopSend,    turbViscBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    3, PETSC_COMM_WORLD, &send_requestT);
    MPI_Irecv(_turbViscBottomRecv, turbViscBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 3, PETSC_COMM_WORLD, &recv_requestBt);

    // Bottom <-- + <-- Top
    MPI_Isend(_turbViscBottomSend, turbViscBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 4, PETSC_COMM_WORLD, &send_requestBt);
    MPI_Irecv(_turbViscTopRecv,    turbViscBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    4, PETSC_COMM_WORLD, &recv_requestT);

    // Wait for all the requests of the current rank to finish before reading
    MPI_Wait(&recv_requestBt, &comm_status);
    MPI_Wait(&recv_requestT,  &comm_status);

    // Read the turbulent viscosity receive buffers T/B
    _parallelBoundaryTurbViscReadIterator->iterate(1);
    // ------------------------------------------------------------------------



    // ------------------------------------------------------------------------
    // Communicate in the Z-dimension  (z increases from front to back):
    if (_parameters.geometry.dim == 3) {

        // Fill the turbulent viscosity send buffers F/B
        _parallelBoundaryTurbViscFillIterator->iterate(2);

        // Front --> + --> Back
        MPI_Isend(_turbViscBackSend,  turbViscBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  5, PETSC_COMM_WORLD, &send_requestBk);
        MPI_Irecv(_turbViscFrontRecv, turbViscBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 5, PETSC_COMM_WORLD, &recv_requestF);

        // Front <-- + <-- Back
        MPI_Isend(_turbViscFrontSend, turbViscBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 6, PETSC_COMM_WORLD, &send_requestF);
        MPI_Irecv(_turbViscBackRecv,  turbViscBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  6, PETSC_COMM_WORLD, &recv_requestBk);

        MPI_Wait(&recv_requestF,  &comm_status);
        MPI_Wait(&recv_requestBk, &comm_status);

        // Read the turbulent viscosity receive buffers F/B
        _parallelBoundaryTurbViscReadIterator->iterate(2);
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
PetscTurbulentParallelManager::~PetscTurbulentParallelManager()  {
    // Free the turbulent viscosity buffers
    free(_turbViscLeftSend);
    free(_turbViscRightSend);
    free(_turbViscBottomSend);
    free(_turbViscTopSend);

    free(_turbViscLeftRecv);
    free(_turbViscRightRecv);
    free(_turbViscBottomRecv);
    free(_turbViscTopRecv);

    if (_parameters.geometry.dim == 3) {
        free(_turbViscFrontSend);
        free(_turbViscBackSend);
        free(_turbViscFrontRecv);
        free(_turbViscBackRecv);
    }
}
