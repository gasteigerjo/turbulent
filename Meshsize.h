#ifndef _MESHSIZE_H_
#define _MESHSIZE_H_

#include "Definitions.h"
#include <cmath>

// forward declaration of Parameters
class Parameters;

enum MeshsizeType{Uniform=0,TanhStretching=1,BfsStretching=2};

/** defines the local mesh size.
 *  @author Philipp Neumann
 */
class Meshsize {
  public:
    virtual ~Meshsize(){}
    // returns the meshsize of cell i,j or i,j,k, respectively
    virtual FLOAT getDx(int i, int j) const = 0;
    virtual FLOAT getDy(int i, int j) const = 0;

    virtual FLOAT getDx(int i, int j, int k) const = 0;
    virtual FLOAT getDy(int i, int j, int k) const = 0;
    virtual FLOAT getDz(int i, int j, int k) const = 0;

    // returns the global geometric position in x-,y-,z-direction
    // of the lower/left/front corner of the local cell at (i,j,k)
    virtual FLOAT getPosX(int i, int j, int k) const = 0;
    virtual FLOAT getPosY(int i, int j, int k) const = 0;
    virtual FLOAT getPosZ(int i, int j, int k) const = 0;

    virtual FLOAT getPosX(int i, int j) const = 0;
    virtual FLOAT getPosY(int i, int j) const = 0;

    // returns the min. meshsize used in this simulation
    // -> required for adaptive time stepping
    virtual FLOAT getDxMin() const = 0;
    virtual FLOAT getDyMin() const = 0;
    virtual FLOAT getDzMin() const = 0;
};


/** implements a uniform, equidistant grid spacing */
class UniformMeshsize: public Meshsize {
  public:
    UniformMeshsize(const Parameters &parameters);
    virtual ~UniformMeshsize();

    virtual FLOAT getDx(int i, int j) const { return _dx; }
    virtual FLOAT getDy(int i, int j) const { return _dy; }

    virtual FLOAT getDx(int i, int j, int k) const { return _dx; }
    virtual FLOAT getDy(int i, int j, int k) const { return _dy; }
    virtual FLOAT getDz(int i, int j, int k) const { return _dz; }

    virtual FLOAT getPosX(int i,int j, int k) const { return _dx*(_firstCornerX-2+i); }
    virtual FLOAT getPosY(int i,int j, int k) const { return _dy*(_firstCornerY-2+j); }
    virtual FLOAT getPosZ(int i,int j, int k) const { return _dz*(_firstCornerZ-2+k); }
    virtual FLOAT getPosX(int i,int j) const { return getPosX(i,j,0); }
    virtual FLOAT getPosY(int i,int j) const { return getPosY(i,j,0); }


    virtual FLOAT getDxMin() const {return _dx;}
    virtual FLOAT getDyMin() const {return _dy;}
    virtual FLOAT getDzMin() const {return _dz;}

  private:
    const FLOAT _dx;
    const FLOAT _dy;
    const FLOAT _dz;
    const int _firstCornerX;
    const int _firstCornerY;
    const int _firstCornerZ;
};


class TanhMeshStretching: public Meshsize {
  public:
    TanhMeshStretching(const Parameters & parameters,bool stretchX, bool stretchY, bool stretchZ);
    virtual ~TanhMeshStretching();

    virtual FLOAT getDx(int i,int j,int k) const {
      if (_stretchX){
        return _coordinatesX[i+1] - _coordinatesX[i];
      } else { return _uniformMeshsize.getDx(i,j,k); }
    }
    virtual FLOAT getDy(int i,int j,int k) const {
      if (_stretchY){
        return _coordinatesY[j+1] - _coordinatesY[j];
      } else { return _uniformMeshsize.getDy(i,j,k); }
    }
    virtual FLOAT getDz(int i,int j,int k) const {
      if (_stretchZ){
        return _coordinatesZ[k+1] - _coordinatesZ[k];
      } else { return _uniformMeshsize.getDz(i,j,k); }
    }
    virtual FLOAT getDx(int i,int j) const { return getDx(i,j,0); }
    virtual FLOAT getDy(int i,int j) const { return getDy(i,j,0); }

    virtual FLOAT getPosX(int i,int j, int k) const {
      if (_stretchX){
        return _coordinatesX[i];
      } else { return _uniformMeshsize.getPosX(i,j,k); }
    }
    virtual FLOAT getPosY(int i,int j, int k) const {
      if (_stretchY){
        return _coordinatesY[j];
      } else { return _uniformMeshsize.getPosY(i,j,k); }
    }
    virtual FLOAT getPosZ(int i,int j, int k) const {
      if (_stretchZ){
        return _coordinatesZ[k];
      } else { return _uniformMeshsize.getPosZ(i,j,k); }
    }
    virtual FLOAT getPosX(int i,int j) const { return getPosX(i,j,0); }
    virtual FLOAT getPosY(int i,int j) const { return getPosY(i,j,0); }


    virtual FLOAT getDxMin() const { return _dxMin; }
    virtual FLOAT getDyMin() const { return _dyMin; }
    virtual FLOAT getDzMin() const { return _dzMin; }


    virtual void precomputeCoordinates();

  protected:
    // computes the coordinate of the lower/left/front corner of the 1D-cell at index i w.r.t. having "size" cells along
    // an interval of length "length". We refer to local indexing, so "firstCorner" denotes the first non-ghost cell index
    // of this process. We use a stretched mesh for all nodes inside the comput. bounding box, and a regular mesh outside this box,
    // using the meshsize of the next inner cell.
    virtual FLOAT computeCoordinate(int i, int firstCorner,int size, FLOAT length, FLOAT deltaS, FLOAT dxMin) const;
    virtual FLOAT computeCoordinateX(int i) const;
    virtual FLOAT computeCoordinateY(int i) const;
    virtual FLOAT computeCoordinateZ(int i) const;

    const Parameters & _parameters;
    const UniformMeshsize _uniformMeshsize;
    const FLOAT _lengthX;
    const FLOAT _lengthY;
    const FLOAT _lengthZ;
    const int _sizeX;
    const int _sizeY;
    const int _sizeZ;
    const int _firstCornerX;
    const int _firstCornerY;
    const int _firstCornerZ;
    const bool _stretchX;
    const bool _stretchY;
    const bool _stretchZ;
    const FLOAT _dxMin;
    const FLOAT _dyMin;
    const FLOAT _dzMin;
    FLOAT * _coordinatesX;
    FLOAT * _coordinatesY;
    FLOAT * _coordinatesZ;
};


class BfsMeshStretching: public TanhMeshStretching {
  public:
    BfsMeshStretching(const Parameters & parameters);
    virtual ~BfsMeshStretching();

    virtual FLOAT getDxMin() const { return std::min(_dxMinBefore, _dxMinAfter); }
    virtual FLOAT getDyMin() const { return std::min(_dyMinBelow, _dyMinAbove); }

    virtual void precomputeCoordinates();

  protected:
    // computes the coordinate of the lower/left/front corner of the 1D-cell at index i w.r.t. having "size" cells along
    // an interval of length "length". We refer to local indexing, so "firstCorner" denotes the first non-ghost cell index
    // of this process. We use a stretched mesh for all nodes inside the comput. bounding box, and a regular mesh outside this box,
    // using the meshsize of the next inner cell.
    virtual FLOAT computeCoordinateX(int i) const;
    virtual FLOAT computeCoordinateY(int i) const;

    const FLOAT _stepX;
    int _sizeXBeforeStep;
    int _sizeXAfterStep;
    FLOAT _dxMinBefore;
    FLOAT _dxMinAfter;
    int _sizeYBelowStep;
    int _sizeYAboveStep;
    FLOAT _dyMinBelow;
    FLOAT _dyMinAbove;
};
#endif // _MESHSIZE_H_
