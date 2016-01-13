#ifndef _CHECKPOINT_H_
#define _CHECKPOINT_H_

#include "Definitions.h"
#include "Parameters.h"
#include "FlowField.h"

class Checkpoint {

    private:

        FlowField & _flowField;
        const Parameters & _parameters;

        MPI_Datatype _filetype;

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        Checkpoint ( FlowField & flowField, const Parameters & parameters );

        ~Checkpoint();

        /** Creates checkpoint files
         * @param timeStep the current timestep
         * @param time the current time
         */
        void create ( int timeStep, FLOAT time );

};

#endif
