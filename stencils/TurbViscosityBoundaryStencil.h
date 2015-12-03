#ifndef _TURB_VISCOSITY_BOUNDARY_STENCIL_H_
#define _TURB_VISCOSITY_BOUNDARY_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../TurbulentFlowField.h"

/**
 * A stencil to set the input velocity in channel flows. Only implements the applyLeftWall(...) methods.
 */
class BFInputTurbViscosityStencil : public BoundaryStencil<TurbulentFlowField> {

    public:
        BFInputTurbViscosityStencil (const Parameters & parameters);

        void applyLeftWall   ( TurbulentFlowField & turbFlowField, int i, int j );
        void applyRightWall  ( TurbulentFlowField & turbFlowField, int i, int j );
        void applyBottomWall ( TurbulentFlowField & turbFlowField, int i, int j );
        void applyTopWall    ( TurbulentFlowField & turbFlowField, int i, int j );

        void applyLeftWall   ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyRightWall  ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyBottomWall ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyTopWall    ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyFrontWall  ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyBackWall   ( TurbulentFlowField & turbFlowField, int i, int j, int k );

    private:
        FLOAT _stepSize; //! fixes the size of the step. If zero, is just channel flow
};

class NeumannTurbViscosityBoundaryStencil: public BoundaryStencil<TurbulentFlowField> {

    public:

        /** Constructor
         * @param parameters Parameters of the simulation
         */
        NeumannTurbViscosityBoundaryStencil(const Parameters & parameters);

        //@brief Functions for the 2D problem. Coordinates entered in alphabetical order.
        //@{
        void applyLeftWall   ( TurbulentFlowField & turbFlowField, int i, int j );
        void applyRightWall  ( TurbulentFlowField & turbFlowField, int i, int j );
        void applyBottomWall ( TurbulentFlowField & turbFlowField, int i, int j );
        void applyTopWall    ( TurbulentFlowField & turbFlowField, int i, int j );
        //@}

        //@brief Functions for the 3D problem. Coordinates entered in alphabetical order.
        //@{
        void applyLeftWall   ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyRightWall  ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyBottomWall ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyTopWall    ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyFrontWall  ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyBackWall   ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        //@}
};

class MovingWallTurbViscosityStencil: public BoundaryStencil<TurbulentFlowField> {

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        MovingWallTurbViscosityStencil ( const Parameters & parameters );

        //@brief Functions for the 2D problem. Coordinates entered in alphabetical order.
        //@{
        void applyLeftWall   ( TurbulentFlowField & turbFlowField, int i, int j );
        void applyRightWall  ( TurbulentFlowField & turbFlowField, int i, int j );
        void applyBottomWall ( TurbulentFlowField & turbFlowField, int i, int j );
        void applyTopWall    ( TurbulentFlowField & turbFlowField, int i, int j );
        //@}

        //@brief Functions for the 3D problem. Coordinates entered in alphabetical order.
        //@{
        void applyLeftWall   ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyRightWall  ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyBottomWall ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyTopWall    ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyFrontWall  ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        void applyBackWall   ( TurbulentFlowField & turbFlowField, int i, int j, int k );
        //@}

};

#endif
