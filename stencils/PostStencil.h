#ifndef _POST_STENCIL_H_
#define _POST_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include <string>
#include <ostream>
#include <sstream>

/** WS1: Stencil for writting VTK files
 *
 * When iterated with, creates a VTK file.
 */
template<class FlowField>
class PostStencil : public FieldStencil<FlowField> {
    private:

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        PostStencil ( const Parameters & parameters )
          : FieldStencil<FlowField>(parameters) {}

        /** Operations to perform before the apply methods
         *
         * @param flowField State of the flow field
         */
        virtual void preapply ( FlowField & flowField ) {} // not pure since it is optional

        /** Writes the information to the file
         * @param stream Stream to be written to
         */
        virtual void writeVtk ( std::ostream & stream ) = 0;

};


class BasicPostStencil : public PostStencil<FlowField> {
    private:
        std::ostringstream _ssPressure;
        std::ostringstream _ssVelocity;
        // std::ostringstream _ssFlags;

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        BasicPostStencil ( const Parameters & parameters );

        /** Operations to perform before the apply methods
         *
         * @param flowField State of the flow field
         */
        void preapply ( FlowField & flowField ) { }

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
         * @param stream Stream to be written to
         */
        void writeVtk ( std::ostream & stream );

};

class WallPostStencil : public PostStencil<FlowField> {
    private:
        std::ostringstream _ssTauw;

        const FLOAT _sizeX, _sizeY, _sizeZ;
        const FLOAT _stepX, _stepY;

        int _jOnStep;
        int _iBehindStep;

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        WallPostStencil ( const Parameters & parameters );

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
         * @param stream Stream to be written to
         */
        void writeVtk ( std::ostream & stream );

};

#endif
