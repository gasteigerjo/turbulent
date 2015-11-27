#ifndef _TURB_VISCOSITY_STENCIL_H_
#define _TURB_VISCOSITY_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../TurbulentFlowField.h"

/** Stencil to compute the turbulent viscosity.
 */
class TurbViscosityStencil : public FieldStencil<TurbulentFlowField> {

    private:

        // A local velocity variable that will be used to approximate derivatives. Size matches 3D
        // case, but can be used for 2D as well.
        FLOAT _localVelocity [ 27 * 3 ];
        // local meshsize
        FLOAT _localMeshsize [ 27 * 3 ];
        
    public:

        /** Constructor
         * @param parameters Parameters of the problem
         */
        TurbViscosityStencil(const Parameters & parameters);

        /** Apply the stencil in 2D
         * @param turbFlowField Turbulent flow field information
         * @param i Position in the X direction
         * @param j Position in the Y direction
         */
        void apply ( TurbulentFlowField & turbFlowField, int i, int j );

        /** Apply the stencil in 3D
         * @param turbFlowField Turbulent flow field information
         * @param i Position in the X direction
         * @param j Position in the Y direction
         * @param k Position in the Z direction
         */
        void apply ( TurbulentFlowField & turbFlowField, int i, int j, int k );
};

#endif
