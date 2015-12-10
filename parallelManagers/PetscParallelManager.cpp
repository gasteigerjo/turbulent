#include "PetscParallelManager.h"

// TODO: Open issues and notes
// - Are the offsets set correctly? We really need to check this!
// - We needed to initialize the stencils and iterators before the body of the constructor. Do some discussion.
// - Is there any problem with the calloc-free method, instead of the new-delete?
// - Iterate the fill and read methods only on the specific boundaries needed.

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

        int lowOffset = 0;//1;//2;
        int highOffset = 0;//-1;

        // int lowOffset = 2;
        // int highOffset = -1;

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
    // int count;

    // ------------------------------------------------------------------------
    // Communicate in the X-dimension:
    // Fill the pressure send buffers L/R
    _parallelBoundaryPressureFillIterator->iterate(0);

    // Left --> + --> Right
    // printf("Communicating L-->R... (rank %d)\n", _parameters.parallel.rank);
    MPI_Sendrecv(_pressuresRightSend, 2*presBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 1,
                 _pressuresLeftRecv,  2*presBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  1,
                 PETSC_COMM_WORLD, &comm_status);
    // MPI_Get_count(&comm_status,MY_MPI_FLOAT,&count);
    // printf("Communicated L-->R (rank %d), count %d from %d\n", _parameters.parallel.rank, count, comm_status.MPI_SOURCE);

    // Left <-- + <-- Right
    // printf("Communicating L<--R... (rank %d)\n", _parameters.parallel.rank);
    MPI_Sendrecv(_pressuresLeftSend,  presBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  1,
                 _pressuresRightRecv, presBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 1,
                 PETSC_COMM_WORLD, &comm_status);
    // MPI_Get_count(&comm_status,MY_MPI_FLOAT,&count);
    // printf("Communicated L<--R (rank %d), count %d from %d\n", _parameters.parallel.rank, count, comm_status.MPI_SOURCE);

    // Wait for the current dimension to complete before continuing to the next
    MPI_Barrier(PETSC_COMM_WORLD);

    // Read the pressure receive buffers L/R
    _parallelBoundaryPressureReadIterator->iterate(0);
    // ------------------------------------------------------------------------



    // ------------------------------------------------------------------------
    // Communicate in the Y-dimension:
    // Fill the pressure send buffers T/B
    _parallelBoundaryPressureFillIterator->iterate(1);

    // Bottom --> + --> Top
    // printf("Communicating B-->T... (rank %d)\n", _parameters.parallel.rank);
    MPI_Sendrecv(_pressuresTopSend,    2*presBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    1,
                 _pressuresBottomRecv, 2*presBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 1,
                 PETSC_COMM_WORLD, &comm_status);
    // MPI_Get_count(&comm_status,MY_MPI_FLOAT,&count);
    // printf("Communicated B-->T (rank %d), count %d from %d\n", _parameters.parallel.rank, count, comm_status.MPI_SOURCE);

    // Bottom <-- + <-- Top
    // printf("Communicating T-->B... (rank %d)\n", _parameters.parallel.rank);
    MPI_Sendrecv(_pressuresBottomSend, presBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 1,
                 _pressuresTopRecv,    presBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    1,
                 PETSC_COMM_WORLD, &comm_status);
    // MPI_Get_count(&comm_status,MY_MPI_FLOAT,&count);
    // printf("Communicated B<--T (rank %d), count %d from %d\n", _parameters.parallel.rank, count, comm_status.MPI_SOURCE);

    // Wait for the current dimension to complete before continuing to the next
    MPI_Barrier(PETSC_COMM_WORLD);

    // Read the pressure receive buffers T/B
    _parallelBoundaryPressureReadIterator->iterate(1);
    // ------------------------------------------------------------------------



    // ------------------------------------------------------------------------
    // Communicate in the Z-dimension  (z increases from front to back):
    if (_parameters.geometry.dim == 3) {

        // Fill the pressure send buffers F/B
        _parallelBoundaryPressureFillIterator->iterate(2);

        // Front --> + --> Back
        // printf("Communicating F-->B (rank %d)\n", _parameters.parallel.rank);
        MPI_Sendrecv(_pressuresBackSend,  2*presBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  1,
                     _pressuresFrontRecv, 2*presBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 1,
                     PETSC_COMM_WORLD, &comm_status);
        // MPI_Get_count(&comm_status,MY_MPI_FLOAT,&count);
        // printf("Communicated F-->B (rank %d), count %d from %d\n", _parameters.parallel.rank, count, comm_status.MPI_SOURCE);

        // Front <-- + <-- Back
        // printf("Communicating F<--B (rank %d)\n", _parameters.parallel.rank);
        MPI_Sendrecv(_pressuresFrontSend, presBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 1,
                     _pressuresBackRecv,  presBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  1,
                     PETSC_COMM_WORLD, &comm_status);
        // MPI_Get_count(&comm_status,MY_MPI_FLOAT,&count);
        // printf("Communicated F<--B (rank %d), count %d from %d\n", _parameters.parallel.rank, count, comm_status.MPI_SOURCE);

        // Wait for the current dimension to complete before continuing to the next
        MPI_Barrier(PETSC_COMM_WORLD);

        // Read the pressure receive buffers F/B
        _parallelBoundaryPressureReadIterator->iterate(2);
    }
    // ------------------------------------------------------------------------

}

void PetscParallelManager::communicateVelocities(){

    MPI_Status comm_status;

    // ------------------------------------------------------------------------
    // Communicate in the X-dimension:
    // Fill the velocity send buffers L/R
    _parallelBoundaryVelocityFillIterator->iterate(0);

    // Left --> + --> Right
    MPI_Sendrecv(_velocitiesRightSend, 2*velBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 1,
                 _velocitiesLeftRecv,  2*velBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  1,
                 PETSC_COMM_WORLD, &comm_status);

    // Left <-- + <-- Right
    MPI_Sendrecv(_velocitiesLeftSend,  velBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.leftNb,  1,
                 _velocitiesRightRecv, velBufSizeLR, MY_MPI_FLOAT, _parameters.parallel.rightNb, 1,
                 PETSC_COMM_WORLD, &comm_status);

    // Wait for the current dimension to complete before continuing to the next
    MPI_Barrier(PETSC_COMM_WORLD);

    // Read the velocity receive buffers L/R
    _parallelBoundaryVelocityReadIterator->iterate(0);
    // ------------------------------------------------------------------------



    // ------------------------------------------------------------------------
    // Communicate in the Y-dimension:
    // Fill the velocity send buffers T/B
    _parallelBoundaryVelocityFillIterator->iterate(1);

    // Bottom --> + --> Top
    MPI_Sendrecv(_velocitiesTopSend,    2*velBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    1,
                 _velocitiesBottomRecv, 2*velBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 1,
                 PETSC_COMM_WORLD, &comm_status);

    // Bottom <-- + <-- Top
    MPI_Sendrecv(_velocitiesBottomSend, velBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.bottomNb, 1,
                 _velocitiesTopRecv,    velBufSizeTB, MY_MPI_FLOAT, _parameters.parallel.topNb,    1,
                 PETSC_COMM_WORLD, &comm_status);

    // Wait for the current dimension to complete before continuing to the next
    MPI_Barrier(PETSC_COMM_WORLD);

    // Read the velocity receive buffers T/B
    _parallelBoundaryVelocityReadIterator->iterate(1);
    // ------------------------------------------------------------------------



    // ------------------------------------------------------------------------
    // Communicate in the Z-dimension (z increases from front to back):
    if (_parameters.geometry.dim == 3) {

        // Fill the pressure send buffers F/B
        _parallelBoundaryVelocityFillIterator->iterate(2);

        // Front --> + --> Back
        MPI_Sendrecv(_velocitiesBackSend,  2*velBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  1,
                     _velocitiesFrontRecv, 2*velBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 1,
                     PETSC_COMM_WORLD, &comm_status);

        // Front <-- + <-- Back
        MPI_Sendrecv(_velocitiesFrontSend, velBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.frontNb, 1,
                     _velocitiesBackRecv,  velBufSizeFB, MY_MPI_FLOAT, _parameters.parallel.backNb,  1,
                     PETSC_COMM_WORLD, &comm_status);

        // Wait for the current dimension to complete
        MPI_Barrier(PETSC_COMM_WORLD);

        // Read the velocity receive buffers F/B
        _parallelBoundaryVelocityReadIterator->iterate(2);
    }
    // ------------------------------------------------------------------------

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
