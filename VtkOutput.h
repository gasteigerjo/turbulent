#ifndef _VTK_OUTPUT_H_
#define _VTK_OUTPUT_H_

#include "Definitions.h"
#include "Parameters.h"
#include "Stencil.h"
#include "FlowField.h"
#include "Iterators.h"
#include "TurbulentFlowField.h"
#include "stencils/PostStencil.h"

/** WS1: Stencil for writing VTK files
 *
 */
class VtkOutput : private FieldStencil<FlowField> {
    private:
        bool _turbulent;

        FlowField & _flowField;
        TurbulentFlowField * _turbFlowField;

        FieldIterator<FlowField> _postIterator;

        int _nPostStencils;
        int _nTurbPostStencils;

        PostStencil<FlowField> ** _postStencils;
        PostStencil<TurbulentFlowField> ** _turbPostStencils;

        void apply ( FlowField & flowField, int i, int j );
        void apply ( FlowField & flowField, int i, int j, int k);

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        VtkOutput ( FlowField & flowField, const Parameters & parameters );

        ~VtkOutput();


        /** Writes the information to the file
         * @param flowField Flow field to be written
         */
        void write ( int timeStep );

};

#endif
