/*----------------------------------------------------------------------------

  CryDisGen, short for Crystal Dislocation Generator, is a software tool
  designed to create general dislocations in crystal structure.

  It is origionally developed by Jianjun Bian, et. al. from Xi'an University 
  of Architeture and Technology and Xi'an Jiaotong Univerisy. 

  Copyright (2025) Xi'an University of Architeture and Technology 
  and Xi'an Jiaotong Univerisy.

  This software is distributed under the GNU General Public License (GPLv3.0). 
  If you have any questions on the usage, please contact the authors at: 
  CryDisGen@gmail.com
  
  ----------------------------------------------------------------------------*/

#ifndef VECT_H
#define VECT_H

typedef struct {
  double x, y, z;
}Vect;

typedef struct {
  int x, y, z;
}VectI;

//direction in lattice
typedef struct {
  double il, im, in;          //Miller index
  double h, k, i, l;          //Miller-Bravais index
  Vect v;                     //vector in Cartesian coodinate system
}Direct;

typedef struct {
  int  id;
  int  type;
  Vect coord;
  int  stat;
}Atom;

typedef struct {
  // v1,v2 are bottom vertices 
  Vect v0, v1, v2; 
}Triangle;

typedef struct {
  Vect     v0, v1;
  Vect     vm0, vm1;
  Triangle t0, t1;
}Dihedral;

typedef struct {
  Vect v0, v1, v2, v3;
  //v0,v1,v3 compose the bottom face
  Dihedral dh0, dh1, dh2;
}Tetrahedron;

typedef struct {
  double xlo, xhi, xprd;
  double ylo, yhi, yprd;
  double zlo, zhi, zprd;
}Region;

typedef struct {
  int i, j;
}Nebr;

// label of the input parameter 
// 0: input exists; 
// 1: no input 
typedef struct {
  int region;
  int x, y, z;
  int lat, pbc;
}RegionParaLabel;

typedef struct {
  int bv, bs;
  int dnorm;
  int dc, dr;
  int lx, ly;
  int ndisl;
}LoopParaLabel;

typedef struct {
  int tc, tr;
  int tnorm;
  int lxt, lyt;
}TetraParaLabel;

typedef struct {
  int bvh, bsh;
  int f;
}HelixParaLabel;

#endif

