#ifndef _DIST_NEAREST_WALL_STENCIL_H_
#define _DIST_NEAREST_WALL_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../TurbulentFlowField.h"

/** Stencil to compute the distance to the nearest wall.
 */
class DistNearestWallStencil : public FieldStencil<TurbulentFlowField> {

    private:
        const FLOAT sizeY;
        const FLOAT sizeZ;
        const FLOAT stepX;
        const FLOAT stepY;

    public:

        /** Constructor
         * @param parameters Parameters of the problem
         */
        DistNearestWallStencil(const Parameters & parameters);

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
