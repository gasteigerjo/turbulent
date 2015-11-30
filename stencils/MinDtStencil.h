#ifndef _MIN_DT_STENCIL_H_
#define _MIN_DT_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../TurbulentFlowField.h"


/**
 */
class MinDtStencil : public FieldStencil<TurbulentFlowField> {

    private:

        FLOAT _minValue;  //!

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        MinDtStencil (const Parameters & parameters);

        void apply (TurbulentFlowField & turbulentFlowField, int i, int j);
        void apply (TurbulentFlowField & turbulentFlowField, int i, int j, int k);

        void reset ();

        /**
         */
        const FLOAT getMinValue();
};

#endif
