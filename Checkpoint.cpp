#include "Checkpoint.h"

Checkpoint::Checkpoint  ( FlowField & flowField, const Parameters & parameters ) :
_flowField(flowField),
_parameters(parameters)
{
    int sizes[3];
    sizes[0] = _parameters.geometry.sizeX;
    sizes[1] = _parameters.geometry.sizeY;
    sizes[2] = 4*_parameters.geometry.sizeZ;

    // Note: we are using a 3D subarray, and we store both the pressure and the three velocity
    //       components in the same subarray, in the form (i,j,k)...(i,j,k+3) for the cell i,j,k.

    int subsizes[3];
    subsizes[0] = _parameters.parallel.localSize[0];
    subsizes[1] = _parameters.parallel.localSize[1];
    subsizes[2] = 4*_parameters.parallel.localSize[2];

    int starts[3];
    starts[0] = _parameters.parallel.firstCorner[0];
    starts[1] = _parameters.parallel.firstCorner[1];
    starts[2] = 4*_parameters.parallel.firstCorner[2];

    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MY_MPI_FLOAT, &_filetype);
}

Checkpoint::~Checkpoint () {}

void Checkpoint::create ( int timeStep, FLOAT time ) {
    MPI_File fh;
    MPI_Status status;
    MPI_Offset disp;
    int ierr;

    ierr = MPI_File_open(PETSC_COMM_WORLD, "my_dummy_file", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    if (ierr != MPI_SUCCESS) {
        handleError(1, "Cannot open/create checkpoint file.");
    }

    if (_parameters.parallel.rank == 0) {

        // ierr = MPI_File_write(fh, &timeStep, 1, MPI_INT, &status);
        //
        // if (ierr != MPI_SUCCESS) {
        //     handleError(1, "Cannot write to the checkpoint file.");
        // }

        ierr = MPI_File_write(fh, &time, 1, MY_MPI_FLOAT, &status);
        // ierr = MPI_File_write(fh, &_flowField.getPressure().getScalar(17,5,5), 1, MY_MPI_FLOAT, &status);

        if (ierr != MPI_SUCCESS) {
            handleError(1, "Cannot write to the checkpoint file.");
        }

        MPI_Offset file_offset;

        // MPI_File_get_position(fh, &file_offset);
        // MPI_File_get_byte_offset(fh, file_offset, &disp);

    }

    // DEBUG: Pray that MPI_INT is legit here instead of MPI_Offset.
    // MPI_Bcast(&disp, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    disp = 8;
    std::cout << "disp = " << disp << std::endl;

    // Note: we assume that we are only using the same machine always, for performance reasons.
    MPI_File_set_view(fh, disp, MY_MPI_FLOAT, _filetype, "internal", MPI_INFO_NULL);

    FLOAT subarray[_parameters.parallel.localSize[0]][_parameters.parallel.localSize[1]][4*_parameters.parallel.localSize[2]];

    // DEBUG: Is the index order optimal??
    for (int i=0; i < _parameters.parallel.localSize[0]; i++) {
        for (int j=0; j < _parameters.parallel.localSize[1]; j++) {
            for (int k=0; k < _parameters.parallel.localSize[2]; k++) {
                // std::cout << _flowField.getPressure().getScalar(i+2,j+2,k+2) << std::endl;
                // printf("%f\n", _flowField.getPressure().getScalar(i+2,j+2,k+2));
                subarray[i][j][4*k] = _flowField.getPressure().getScalar(i+2,j+2,k+2);
                subarray[i][j][4*k+1] = _flowField.getVelocity().getVector(i+2,j+2,k+2)[0];
                subarray[i][j][4*k+2] = _flowField.getVelocity().getVector(i+2,j+2,k+2)[1];
                subarray[i][j][4*k+3] = _flowField.getVelocity().getVector(i+2,j+2,k+2)[2];
            }
        }
    }

    int data_count = 4 * _parameters.parallel.localSize[0] * _parameters.parallel.localSize[1] * _parameters.parallel.localSize[2];

    ierr = MPI_File_write(fh, subarray, data_count, MY_MPI_FLOAT, &status);

    // int mycount;
    // MPI_Get_count(&status, MY_MPI_FLOAT, &mycount);
    // std::cout << "count = " << mycount << std::endl;
    // std::cout << "undefined = " << MPI_UNDEFINED << std::endl;

    if (ierr != MPI_SUCCESS) {
        handleError(1, "Cannot write to the checkpoint file.");
    }

    MPI_File_close(&fh);
}
