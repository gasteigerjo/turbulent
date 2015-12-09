#ifndef _VTK_STENCIL_H_
#define _VTK_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include <string>
#include <fstream>
#include <sstream>

/** WS1: Stencil for writting VTK files
 *
 * When iterated with, creates a VTK file.
 */
class VTKStencil : public FieldStencil<FlowField> {
    private:
        std::stringstream _ssPoints, _ssPressure, _ssFlags, _ssVelocity, _ssTurbViscosity, _ssWallDistance;
        bool _turbulent;
        bool _includeGhostCells;

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        VTKStencil ( const Parameters & parameters );

        /** Constructor
         *
         * @param parameters Parameters of the problem
         * @param includeGhostCells Switch to include the ghost layers in the vtk
         */
        VTKStencil ( const Parameters & parameters, bool includeGhostCells );

        /** 2D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         */
        void apply ( FlowField & flowField, int i, int j );

        /** 3D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         * @param k Position in the z direction
         */
        void apply ( FlowField & flowField, int i, int j, int k );

        /** Writes the information to the file
         * @param flowField Flow field to be written
         */
        void write ( FlowField & flowField, int timeStep );

};

#endif
