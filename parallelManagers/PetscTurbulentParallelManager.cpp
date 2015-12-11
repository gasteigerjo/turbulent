#include "PetscTurbulentParallelManager.h"

// TODO: Open issues and notes
// - Are the offsets set correctly? We really need to check this!
// - We needed to initialize the stencils and iterators before the body of the constructor. Do some discussion.
// - Is there any problem with the calloc-free method, instead of the new-delete?
// - Iterate the fill and read methods only on the specific boundaries needed.

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
    // int count;

    // ------------------------------------------------------------------------
    // Communicate in the X-dimension:
    // Fill the turbulent viscosity send buffers L/R
    _parallelBoundaryTurbViscFillIterator->iterate(0);

    // Left --> + --> Right
    // printf("Communicating L-->R... (rank %d)\n", _parameters.parallel.rank);
    MPI_Sendrecv(_turbViscRightSend, turbViscBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 1,
                 _turbViscLeftRecv,  turbViscBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  1,
                 PETSC_COMM_WORLD, &comm_status);
    // MPI_Get_count(&comm_status,MY_MPI_FLOAT,&count);
    // printf("Communicated L-->R (rank %d), count %d from %d\n", _parameters.parallel.rank, count, comm_status.MPI_SOURCE);

    // Left <-- + <-- Right
    // printf("Communicating L<--R... (rank %d)\n", _parameters.parallel.rank);
    MPI_Sendrecv(_turbViscLeftSend,  turbViscBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  1,
                 _turbViscRightRecv, turbViscBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 1,
                 PETSC_COMM_WORLD, &comm_status);
    // MPI_Get_count(&comm_status,MY_MPI_FLOAT,&count);
    // printf("Communicated L<--R (rank %d), count %d from %d\n", _parameters.parallel.rank, count, comm_status.MPI_SOURCE);

    // Wait for the current dimension to complete before continuing to the next
    MPI_Barrier(PETSC_COMM_WORLD);

    // Read the turbulent viscosity receive buffers L/R
    _parallelBoundaryTurbViscReadIterator->iterate(0);
    // ------------------------------------------------------------------------



    // ------------------------------------------------------------------------
    // Communicate in the Y-dimension:
    // Fill the turbulent viscosity send buffers T/B
    _parallelBoundaryTurbViscFillIterator->iterate(1);

    // Bottom --> + --> Top
    // printf("Communicating B-->T... (rank %d)\n", _parameters.parallel.rank);
    MPI_Sendrecv(_turbViscTopSend,    turbViscBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    1,
                 _turbViscBottomRecv, turbViscBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 1,
                 PETSC_COMM_WORLD, &comm_status);
    // MPI_Get_count(&comm_status,MY_MPI_FLOAT,&count);
    // printf("Communicated B-->T (rank %d), count %d from %d\n", _parameters.parallel.rank, count, comm_status.MPI_SOURCE);

    // Bottom <-- + <-- Top
    // printf("Communicating T-->B... (rank %d)\n", _parameters.parallel.rank);
    MPI_Sendrecv(_turbViscBottomSend, turbViscBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 1,
                 _turbViscTopRecv,    turbViscBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    1,
                 PETSC_COMM_WORLD, &comm_status);
    // MPI_Get_count(&comm_status,MY_MPI_FLOAT,&count);
    // printf("Communicated B<--T (rank %d), count %d from %d\n", _parameters.parallel.rank, count, comm_status.MPI_SOURCE);

    // Wait for the current dimension to complete before continuing to the next
    MPI_Barrier(PETSC_COMM_WORLD);

    // Read the turbulent viscosity receive buffers T/B
    _parallelBoundaryTurbViscReadIterator->iterate(1);
    // ------------------------------------------------------------------------



    // ------------------------------------------------------------------------
    // Communicate in the Z-dimension  (z increases from front to back):
    if (_parameters.geometry.dim == 3) {

        // Fill the turbulent viscosity send buffers F/B
        _parallelBoundaryTurbViscFillIterator->iterate(2);

        // Front --> + --> Back
        // printf("Communicating F-->B (rank %d)\n", _parameters.parallel.rank);
        MPI_Sendrecv(_turbViscBackSend,  turbViscBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  1,
                     _turbViscFrontRecv, turbViscBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 1,
                     PETSC_COMM_WORLD, &comm_status);
        // MPI_Get_count(&comm_status,MY_MPI_FLOAT,&count);
        // printf("Communicated F-->B (rank %d), count %d from %d\n", _parameters.parallel.rank, count, comm_status.MPI_SOURCE);

        // Front <-- + <-- Back
        // printf("Communicating F<--B (rank %d)\n", _parameters.parallel.rank);
        MPI_Sendrecv(_turbViscFrontSend, turbViscBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 1,
                     _turbViscBackRecv,  turbViscBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  1,
                     PETSC_COMM_WORLD, &comm_status);
        // MPI_Get_count(&comm_status,MY_MPI_FLOAT,&count);
        // printf("Communicated F<--B (rank %d), count %d from %d\n", _parameters.parallel.rank, count, comm_status.MPI_SOURCE);

        // Wait for the current dimension to complete before continuing to the next
        MPI_Barrier(PETSC_COMM_WORLD);

        // Read the turbulent viscosity receive buffers F/B
        _parallelBoundaryTurbViscReadIterator->iterate(2);
    }
    // ------------------------------------------------------------------------

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
