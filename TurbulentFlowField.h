#ifndef _TURB_FLOW_FIELD_H_
#define _TURB_FLOW_FIELD_H_

#include "DataStructures.h"
#include "Parameters.h"
#include "FlowField.h"

/** Turbulent flow field
 *
 * Class intended to contain the state of the domain.
 */
class TurbulentFlowField : public FlowField {

    private:

        ScalarField _turbViscosity;   //! Scalar field representing the turbulent viscosity
        ScalarField _distNearestWall; //! Scalar field representing the distance to the nearest wall

    public:
        /** Get turbulent viscosity field
         * @return Reference to turbulent viscosity field
         */
        ScalarField & getTurbViscosity ();

        /** Get field with distance to nearest wall
         * @return Reference to field
         */
        ScalarField & getDistNearestWall ();
};

#endif
