#ifndef _PETSC_PARALLEL_MANAGER_H_
#define _PETSC_PARALLEL_MANAGER_H_

#include "../stencils/PressureBufferFillStencil.h"
#include "../stencils/PressureBufferReadStencil.h"
#include "../stencils/VelocityBufferFillStencil.h"
#include "../stencils/VelocityBufferReadStencil.h"
#include "../Parameters.h"
#include "../FlowField.h"
#include "../Iterators.h"
#include "../Definitions.h"

class PetscParallelManager {

    private:

        const Parameters & _parameters;
        FlowField & _flowField;

        // Size of buffers in each dimension
        // sizes of the pressure buffers
        int presBufSizeLR; // For communicating Left<-->Right.
        int presBufSizeTB; // For communicating Top<-->Bottom.
        int presBufSizeFB; // For communicating Front<-->Back.
        // sizes of the velocity buffers
        int velBufSizeLR; // For communicating Left<-->Right.
        int velBufSizeTB; // For communicating Top<-->Bottom.
        int velBufSizeFB; // For communicating Front<-->Back.

        // Send buffers for pressure
        FLOAT * _pressuresLeftSend;
        FLOAT * _pressuresRightSend;
        FLOAT * _pressuresBottomSend;
        FLOAT * _pressuresTopSend;
        FLOAT * _pressuresFrontSend;
        FLOAT * _pressuresBackSend;

        // Receive buffers for pressure
        FLOAT * _pressuresLeftRecv;
        FLOAT * _pressuresRightRecv;
        FLOAT * _pressuresBottomRecv;
        FLOAT * _pressuresTopRecv;
        FLOAT * _pressuresFrontRecv;
        FLOAT * _pressuresBackRecv;

        // Send buffers for velocities
        FLOAT * _velocitiesLeftSend;
        FLOAT * _velocitiesRightSend;
        FLOAT * _velocitiesBottomSend;
        FLOAT * _velocitiesTopSend;
        FLOAT * _velocitiesFrontSend;
        FLOAT * _velocitiesBackSend;

        // Receive buffers for velocities
        FLOAT * _velocitiesLeftRecv;
        FLOAT * _velocitiesRightRecv;
        FLOAT * _velocitiesBottomRecv;
        FLOAT * _velocitiesTopRecv;
        FLOAT * _velocitiesFrontRecv;
        FLOAT * _velocitiesBackRecv;

        // Stencil and iterator objects
        PressureBufferFillStencil * _pressureBufferFillStencil;
        PressureBufferReadStencil * _pressureBufferReadStencil;
        VelocityBufferFillStencil * _velocityBufferFillStencil;
        VelocityBufferReadStencil * _velocityBufferReadStencil;
        ParallelBoundaryIterator<FlowField> * _parallelBoundaryPressureFillIterator;
        ParallelBoundaryIterator<FlowField> * _parallelBoundaryPressureReadIterator;
        ParallelBoundaryIterator<FlowField> * _parallelBoundaryVelocityFillIterator;
        ParallelBoundaryIterator<FlowField> * _parallelBoundaryVelocityReadIterator;

    public:

        /** Constructor
         * @param parameters Reference to the parameters
         */
        PetscParallelManager(const Parameters & parameters, FlowField & flowField);

        /** Destructor */
        ~PetscParallelManager();

        /** Communicates the pressure buffers between processes.
        */
        void communicatePressure();

        /** Communicates the velocity buffers between processes.
        */
        void communicateVelocities();

};


#endif
