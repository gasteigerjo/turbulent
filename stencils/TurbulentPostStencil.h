#ifndef _TURB_POST_STENCIL_H_
#define _TURB_POST_STENCIL_H_

#include "../DataStructures.h"
#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../TurbulentFlowField.h"
#include "PostStencil.h"
#include <string>
#include <ostream>
#include <sstream>

/** WS1: Stencil for writting VTK files
 *
 * When iterated with, creates a VTK file.
 */
class TurbulentPostStencil : public PostStencil<TurbulentFlowField> {
    private:
        std::ostringstream _ssTurbVisc;
        // std::ostringstream _ssDistWall;

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        TurbulentPostStencil ( const Parameters & parameters )
          : PostStencil<TurbulentFlowField>(parameters) {}

        /** Operations to perform before the apply methods
         *
         * @param flowField State of the flow field
         */
        void preapply ( TurbulentFlowField & flowField ) { }

        /** 2D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         */
        void apply ( TurbulentFlowField & flowField, int i, int j );

        /** 3D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         * @param k Position in the z direction
         */
        void apply ( TurbulentFlowField & flowField, int i, int j, int k );

        /** Writes the information to the file
         * @param stream Stream to be written to
         */
        void writeVtk ( std::ostream & stream );

};

class TurbulentWallPostStencil : public PostStencil<TurbulentFlowField> {
    private:
        std::ostringstream _ssUTau;
        std::ostringstream _ssYPlus;
        std::ostringstream _ssUPlus;

        const FLOAT _sizeX, _sizeY, _sizeZ;
        const FLOAT _stepX, _stepY;

        const bool _possible;

        FLOAT * _uTau_bottom;
        FLOAT * _uTau_top;
        FLOAT * _uTau_front;
        FLOAT * _uTau_back;
        FLOAT * _uTau_step;

        int _jOnStep;
        int _iBehindStep;

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        TurbulentWallPostStencil ( const Parameters & parameters );
          // : PostStencil<TurbulentFlowField>(parameters) {}

        ~TurbulentWallPostStencil();

        /** Operations to perform before the apply methods
         *
         * @param flowField State of the flow field
         */
        void preapply ( TurbulentFlowField & flowField );

        /** 2D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         */
        void apply ( TurbulentFlowField & flowField, int i, int j );

        /** 3D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         * @param k Position in the z direction
         */
        void apply ( TurbulentFlowField & flowField, int i, int j, int k );

        /** Writes the information to the file
         * @param stream Stream to be written to
         */
        void writeVtk ( std::ostream & stream );

};

#endif
