#!/bin/bash

# Filename of the binary file to convert.
FILENAME=$1
# Number of cells written in the file.
CELLS=$2
# Filename of the ASCII file to convert to.
OUTPUT=${FILENAME}.txt
# By default, assume 3D case. Only if $3 is given and is equal to 2, change to the 2D case.
if [ "$3" = "2" ]; then
    DIM=2
else
    DIM=3
fi

# Read and write the first 12 bytes of the binary file (an integer and a float).
hexdump -e '1 4 "Timestep: %d \n" 1 8 "Time    : %f\n\n"' -n 12 ${FILENAME} > ${OUTPUT}
# We have already read 12 bytes.
OFFSET=12

# Write the column titles (2D or 3D case).
if [ "$DIM" = "2" ]; then
    printf "%s\n" " Index    Pressure           Velocity [x,y]          " >> ${OUTPUT}
    printf "%s\n" "-----------------------------------------------------" >> ${OUTPUT}
else
    printf "%s\n" " Index    Pressure               Velocity [x,y,z]          " >> ${OUTPUT}
    printf "%s\n" "-----------------------------------------------------------" >> ${OUTPUT}
fi

if [ "$DIM" = "2" ]; then
    # Read and write the rest of the file, 24 bytes (3 floats) per time.
    for i in `seq 1 ${CELLS}`;
    do 
        printf "[%6d] " ${i} >> ${OUTPUT}
        hexdump -e '1 8 "%10.6f    | " 2 8 "%10.6f " "\n"' -n 24 -s ${OFFSET} ${FILENAME} >> ${OUTPUT}
        let OFFSET+=24
    done    
else
    # Read and write the rest of the file, 32 bytes (4 floats) per time.
    for i in `seq 1 ${CELLS}`;
    do 
        printf "[%6d] " ${i} >> ${OUTPUT}
        hexdump -e '1 8 "%10.6f    | " 3 8 "%10.6f " "\n"' -n 32 -s ${OFFSET} ${FILENAME} >> ${OUTPUT}
        let OFFSET+=32
    done
fi
