#ifndef _TURB_VISCOSITY_BUFFER_FILL_STENCIL_H_
#define _TURB_VISCOSITY_BUFFER_FILL_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../TurbulentFlowField.h"
#include "ScalarBufferFillStencil.h"
#include <string>

/**
 *
 * Fills an array with the turbulent viscosity values at the boundary.
 */
class TurbViscosityBufferFillStencil : public ScalarBufferFillStencil<TurbulentFlowField> {

public:
    /** Constructor for 2D case
     * @param parameters Parameters of the simulation
     * @param turbViscLeft Pointer to an array of length (N+2) for saving the turbulent viscosity at the left boundary
     * @param turbViscRight Pointer to an array of length (N+2) for saving the turbulent viscosity at the right boundary
     * @param turbViscBottom Pointer to an array of length (N+2) for saving the turbulent viscosity at the bottom boundary
     * @param turbViscTop Pointer to an array of length (N+2) for saving the turbulent viscosity at the top boundary
     */
    TurbViscosityBufferFillStencil(const Parameters & parameters,
        FLOAT * turbViscLeft, FLOAT * turbViscRight,
        FLOAT * turbViscBottom, FLOAT * turbViscTop,
        int lowOffset = 0);

    /** Constructor for 3D case
     * @param parameters Parameters of the simulation
     * @param turbViscLeft Pointer to an array of length (N+2)^2 for saving the turbulent viscosity at the left boundary
     * @param turbViscRight Pointer to an array of length (N+2)^2 for saving the turbulent viscosity at the right boundary
     * @param turbViscBottom Pointer to an array of length (N+2)^2 for saving the turbulent viscosity at the bottom boundary
     * @param turbViscTop Pointer to an array of length (N+2)^2 for saving the turbulent viscosity at the top boundary
     * @param turbViscFront Pointer to an array of length (N+2)^2 for saving the turbulent viscosity at the front boundary
     * @param turbViscBack Pointer to an array of length (N+2)^2 for saving the turbulent viscosity at the back boundary
     */
    TurbViscosityBufferFillStencil(const Parameters & parameters,
        FLOAT * turbViscLeft, FLOAT * turbViscRight,
        FLOAT * turbViscBottom, FLOAT * turbViscTop,
        FLOAT * turbViscFront, FLOAT * turbViscBack,
        int lowOffset = 0);

protected:
    FLOAT & getScalar( TurbulentFlowField & turbFlowField, int i, int j );
    FLOAT & getScalar( TurbulentFlowField & turbFlowField, int i, int j, int k );

};

#endif
