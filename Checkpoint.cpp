#include "Checkpoint.h"

Checkpoint::Checkpoint  ( FlowField & flowField, const Parameters & parameters ) :
_flowField(flowField),
_parameters(parameters)
{
    // 2D or 3D case?
    if (_parameters.geometry.dim == 2) {
    // 2D case
        int sizes[2];
        sizes[0] = _parameters.geometry.sizeX;
        sizes[1] = 3*_parameters.geometry.sizeY;

        // Note: we are using a 2D subarray, and we store both the pressure and the three velocity
        //       components in the same subarray, in the form (i,j)...(i,j+2) for each cell i,j.
        int subsizes[2];
        subsizes[0] = _parameters.parallel.localSize[0];
        subsizes[1] = 3*_parameters.parallel.localSize[1];

        int starts[2];
        starts[0] = _parameters.parallel.firstCorner[0];
        starts[1] = 3*_parameters.parallel.firstCorner[1];

        MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MY_MPI_FLOAT, &_filetype);
    } else {
    // 3D case
        int sizes[3];
        sizes[0] = _parameters.geometry.sizeX;
        sizes[1] = _parameters.geometry.sizeY;
        sizes[2] = 4*_parameters.geometry.sizeZ;

        // Note: we are using a 3D subarray, and we store both the pressure and the three velocity
        //       components in the same subarray, in the form (i,j,k)...(i,j,k+3) for each cell i,j,k.
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
    MPI_Type_commit(&_filetype);

    // Create the restart directory if it doesn't exist
    mkdir(_parameters.checkpoint.directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH);

}

Checkpoint::~Checkpoint () {}

void Checkpoint::read ( int& timeStep, FLOAT& time ) {
    MPI_File fh_restart;
    MPI_Status status;
    MPI_Offset disp;
    int ierr;

    // Open the restart file
    char* tmp_filename = new char[strlen(_parameters.restart.filename.c_str())];
    strcpy(tmp_filename, _parameters.restart.filename.c_str());
    ierr = MPI_File_open(PETSC_COMM_WORLD, tmp_filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_restart);
    delete[] tmp_filename;
    
    if (ierr != MPI_SUCCESS) {
        handleError(1, "Cannot open the restart file.");
    }

    int buffer_timeStep;
    // Read the timeStep
    MPI_File_read_at(fh_restart, 0, &buffer_timeStep, 1, MPI_INT, &status);
    if (ierr != MPI_SUCCESS) {
        handleError(1, "Cannot read the timeStep from the restart file.");
    }
    // DEBUG
    printf("====== TIMESTEP: %d ==========\n", buffer_timeStep);
    timeStep = buffer_timeStep;

    FLOAT buffer_time;
    // Read the time
    MPI_File_read_at(fh_restart, sizeof(int), &buffer_time, 1, MY_MPI_FLOAT, &status);
    if (ierr != MPI_SUCCESS) {
        handleError(1, "Cannot read the time from the restart file.");
    }
    // DEBUG
    printf("====== TIME: %f ==========\n", buffer_time);
    time = buffer_time;

    // Displacement of the file view from the begining of the file.
    // The header consists of an int and a FLOAT.
    disp = sizeof(int) + sizeof(FLOAT);

    // Note: we assume that we are only using the same machine always, for performance reasons.
    MPI_File_set_view(fh_restart, disp, MY_MPI_FLOAT, _filetype, "native", MPI_INFO_NULL);

    // 2D or 3D case?
    if (_parameters.geometry.dim == 2) {
    // 2D case
        FLOAT localarray[_parameters.parallel.localSize[0]][3*_parameters.parallel.localSize[1]];

        // Number of elements to read
        int data_count = 3 * _parameters.parallel.localSize[0] * _parameters.parallel.localSize[1];

        // Read the cell data from the file
        ierr = MPI_File_read(fh_restart, localarray, data_count, MY_MPI_FLOAT, &status);
        if (ierr != MPI_SUCCESS) {
            handleError(1, "Cannot read the cell data from the checkpoint file.");
        }

        for (int i=0; i < _parameters.parallel.localSize[0]; i++) {
            for (int j=0; j < _parameters.parallel.localSize[1]; j++) {
                // DEBUG
                if (_parameters.parallel.rank == 0) {
                    printf("%f | %f %f\n", localarray[i][3*j], localarray[i][3*j+1], localarray[i][3*j+2]);
                }
                _flowField.getPressure().getScalar(i+2,j+2) = localarray[i][3*j];
                _flowField.getVelocity().getVector(i+2,j+2)[0] = localarray[i][3*j+1];
                _flowField.getVelocity().getVector(i+2,j+2)[1] = localarray[i][3*j+2];
            }
        }
    } else {
    // 3D case
        FLOAT localarray[_parameters.parallel.localSize[0]][_parameters.parallel.localSize[1]][4*_parameters.parallel.localSize[2]];

        // Number of elements to write
        int data_count = 4 * _parameters.parallel.localSize[0] * _parameters.parallel.localSize[1] * _parameters.parallel.localSize[2];

        // Read the cell data from the file
        ierr = MPI_File_read(fh_restart, localarray, data_count, MY_MPI_FLOAT, &status);
        if (ierr != MPI_SUCCESS) {
            handleError(1, "Cannot read the cell data from the checkpoint file.");
        }

        for (int i=0; i < _parameters.parallel.localSize[0]; i++) {
            for (int j=0; j < _parameters.parallel.localSize[1]; j++) {
                for (int k=0; k < _parameters.parallel.localSize[2]; k++) {
                    // DEBUG
                    if (_parameters.parallel.rank == 0) {
                        printf("%f | %f %f %f\n", localarray[i][j][4*k], localarray[i][j][4*k+1], localarray[i][j][4*k+2], localarray[i][j][4*k+3]);
                    }
                    _flowField.getPressure().getScalar(i+2,j+2,k+2) = localarray[i][j][4*k];
                    _flowField.getVelocity().getVector(i+2,j+2,k+2)[0] = localarray[i][j][4*k+1];
                    _flowField.getVelocity().getVector(i+2,j+2,k+2)[1] = localarray[i][j][4*k+2];
                    _flowField.getVelocity().getVector(i+2,j+2,k+2)[2] = localarray[i][j][4*k+3];
                }
            }
        }
    }

    // Close the file
    ierr = MPI_File_close(&fh_restart);
    if (ierr != MPI_SUCCESS) {
        handleError(1, "Cannot close the restart file.");
    }
}

void Checkpoint::create ( int timeStep, FLOAT time ) {
    MPI_File fh_checkpoint;
    MPI_Status status;
    MPI_Offset disp;
    int ierr;

    // Set the filename of the checkpoint file
    std::ostringstream checkpointFileName;
    checkpointFileName << _parameters.checkpoint.directory << _parameters.checkpoint.prefix
             << "." << std::setfill('0') << std::setw(6) << timeStep;

    // Open the file and assign it to the handler fh_checkpoint.
    char* tmp_filename = new char[strlen(checkpointFileName.str().c_str())];
    strcpy(tmp_filename, checkpointFileName.str().c_str());
    ierr = MPI_File_open(PETSC_COMM_WORLD, tmp_filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_checkpoint);
    delete[] tmp_filename;

    if (ierr != MPI_SUCCESS) {
        handleError(1, "Cannot open/create checkpoint file.");
    }

    // Write the header of the file using Rank0.
    if (_parameters.parallel.rank == 0) {

        // Write the timeStep
        ierr = MPI_File_write(fh_checkpoint, &timeStep, 1, MPI_INT, &status);
        if (ierr != MPI_SUCCESS) {
            handleError(1, "Cannot write the timeStep to the checkpoint file.");
        }

        // Write the time
        ierr = MPI_File_write(fh_checkpoint, &time, 1, MY_MPI_FLOAT, &status);
        if (ierr != MPI_SUCCESS) {
            handleError(1, "Cannot write the time to the checkpoint file.");
        }

    }

    // Displacement of the file view from the begining of the file.
    // The header we wrote consists of an int and a FLOAT.
    disp = sizeof(int) + sizeof(FLOAT);

    // Note: we assume that we are only using the same machine always, for performance reasons.
    MPI_File_set_view(fh_checkpoint, disp, MY_MPI_FLOAT, _filetype, "native", MPI_INFO_NULL);

    // 2D or 3D case?
    if (_parameters.geometry.dim == 2) {
    // 2D case
        FLOAT localarray[_parameters.parallel.localSize[0]][3*_parameters.parallel.localSize[1]];

        // Prepare the cell data for writting to the file.
        // DEBUG: Is the index order optimal??
        for (int i=0; i < _parameters.parallel.localSize[0]; i++) {
            for (int j=0; j < _parameters.parallel.localSize[1]; j++) {
                // DEBUG
                // if (_parameters.parallel.rank == 0) {
                //     printf("%f | %f %f\n", _flowField.getPressure().getScalar(i+2,j+2), _flowField.getVelocity().getVector(i+2,j+2)[0], _flowField.getVelocity().getVector(i+2,j+2)[1]);
                // }
                // Pressure
                localarray[i][3*j] = _flowField.getPressure().getScalar(i+2,j+2);
                // Velocities per x,y
                localarray[i][3*j+1] = _flowField.getVelocity().getVector(i+2,j+2)[0];
                localarray[i][3*j+2] = _flowField.getVelocity().getVector(i+2,j+2)[1];
            }
        }

        // Number of elements to write
        int data_count = 3 * _parameters.parallel.localSize[0] * _parameters.parallel.localSize[1];

        // Write the cell data to the file
        ierr = MPI_File_write(fh_checkpoint, localarray, data_count, MY_MPI_FLOAT, &status);
        if (ierr != MPI_SUCCESS) {
            handleError(1, "Cannot write the cell data to the checkpoint file.");
        }
    } else {
    // 3D case
        FLOAT localarray[_parameters.parallel.localSize[0]][_parameters.parallel.localSize[1]][4*_parameters.parallel.localSize[2]];

        // Prepare the data to write to the file.
        // DEBUG: Is the index order optimal??
        for (int i=0; i < _parameters.parallel.localSize[0]; i++) {
            for (int j=0; j < _parameters.parallel.localSize[1]; j++) {
                for (int k=0; k < _parameters.parallel.localSize[2]; k++) {
                    // DEBUG
                    // if (_parameters.parallel.rank == 0) {
                    //     printf("%f | %f %f %f\n", _flowField.getPressure().getScalar(i+2,j+2,k+2), _flowField.getVelocity().getVector(i+2,j+2,k+2)[0], _flowField.getVelocity().getVector(i+2,j+2,k+2)[1], _flowField.getVelocity().getVector(i+2,j+2,k+2)[2]);
                    // }
                    // Pressure
                    localarray[i][j][4*k] = _flowField.getPressure().getScalar(i+2,j+2,k+2);
                    // Velocities per x,y,z
                    localarray[i][j][4*k+1] = _flowField.getVelocity().getVector(i+2,j+2,k+2)[0];
                    localarray[i][j][4*k+2] = _flowField.getVelocity().getVector(i+2,j+2,k+2)[1];
                    localarray[i][j][4*k+3] = _flowField.getVelocity().getVector(i+2,j+2,k+2)[2];
                }
            }
        }

        // Number of elements to write
        int data_count = 4 * _parameters.parallel.localSize[0] * _parameters.parallel.localSize[1] * _parameters.parallel.localSize[2];

        // Write the cell data to the file
        ierr = MPI_File_write(fh_checkpoint, localarray, data_count, MY_MPI_FLOAT, &status);
        if (ierr != MPI_SUCCESS) {
            handleError(1, "Cannot write the cell data to the checkpoint file.");
        }
    }

    // Close the file
    ierr = MPI_File_close(&fh_checkpoint);
    if (ierr != MPI_SUCCESS) {
        handleError(1, "Cannot close the checkpoint file.");
    }
}

void Checkpoint::cleandir () {

    // Clean the directory from previous restart data
    DIR *checkpointDir = opendir(_parameters.checkpoint.directory.c_str());
    struct dirent *next_file;
    char filepath[256];

    while( (next_file = readdir(checkpointDir)) != NULL ) {
        sprintf(filepath, "%s/%s", _parameters.checkpoint.directory.c_str(), next_file->d_name);
        remove(filepath);
    }
    closedir(checkpointDir);

}
