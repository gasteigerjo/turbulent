#ifndef _PETSC_TURBULENT_PARALLEL_MANAGER_H_
#define _PETSC_TURBULENT_PARALLEL_MANAGER_H_

#include "../stencils/TurbViscosityBufferFillStencil.h"
#include "../stencils/TurbViscosityBufferReadStencil.h"
#include "../Parameters.h"
#include "../TurbulentFlowField.h"
#include "../Iterators.h"
#include "../Definitions.h"
#include "PetscParallelManager.h"

class PetscTurbulentParallelManager {//: public PetscParallelManager {

    private:
        const Parameters & _parameters;
        TurbulentFlowField & _turbFlowField;

        // Size of buffers in each dimension
        // sizes of the turbulent viscosity buffers
        int turbViscBufSizeLR; // For communicating Left<-->Right.
        int turbViscBufSizeTB; // For communicating Top<-->Bottom.
        int turbViscBufSizeFB; // For communicating Front<-->Back.

        // Send buffers for turbulent viscosity
        FLOAT * _turbViscLeftSend;
        FLOAT * _turbViscRightSend;
        FLOAT * _turbViscBottomSend;
        FLOAT * _turbViscTopSend;
        FLOAT * _turbViscFrontSend;
        FLOAT * _turbViscBackSend;

        // Receive buffers for turbulent viscosity
        FLOAT * _turbViscLeftRecv;
        FLOAT * _turbViscRightRecv;
        FLOAT * _turbViscBottomRecv;
        FLOAT * _turbViscTopRecv;
        FLOAT * _turbViscFrontRecv;
        FLOAT * _turbViscBackRecv;

        // Stencil and iterator objects
        TurbViscosityBufferFillStencil * _turbViscBufferFillStencil;
        TurbViscosityBufferReadStencil * _turbViscBufferReadStencil;
        ParallelBoundaryIterator<TurbulentFlowField> * _parallelBoundaryTurbViscFillIterator;
        ParallelBoundaryIterator<TurbulentFlowField> * _parallelBoundaryTurbViscReadIterator;

    public:

        /** Constructor
         * @param parameters Reference to the parameters
         */
        PetscTurbulentParallelManager(const Parameters & parameters, TurbulentFlowField & fturbFowField);

        /** Destructor */
        ~PetscTurbulentParallelManager();

        /** Communicates the turbulent viscosity buffers between processes.
        */
        void communicateTurbViscosity();

};


#endif
