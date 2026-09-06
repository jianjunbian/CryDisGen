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

#ifndef CREATEDISL_H
#define CREATEDISL_H

#include "vect.h"
#include <string.h>

class CreateDisl{
  public:
    CreateDisl(int, char**);
    void   ReadParaFile(char*);
    void   ReadConfig(char*);

    //dislocation loop
    void   DiscretizeLoop(int);
    double AtomSolidAngle(int,Vect);
    void   ApplyDisplacement(int);

		//elastic displacement
		Vect   ElasticDispLoop(Vect,int);
		Vect   ElasticDispTriangle(Vect,Vect,Vect,Vect,Vect);
		Vect   ElasticDispHeli(int,Vect,Vect);
		//to avoid atoms perfectly align with dislocation core
		void   ShiftAtoms(); 

    //tetrahedron stacking fault
    void   DiscretizeLoopTe(int);
    double AtomSolidAngleTe(int,int,Vect);
    void   ApplyDisplacementTe(int,int);

    //helical dislocation
    void   ReadGrid(int);
    double AtomSolidAngleHeli(int,Vect);
    void   ApplyDisplacementHeli(int);

    void   Output(char*);

  private:
    //function to calculate solid angle
    double SolidAngleTriangle(Triangle, Vect);
    double SolidAngleTetrahedron(Tetrahedron);
		//calculate elastic displacement at a field point 
		//with respect to a straight dislocation segment
		Vect   DispSegment(Vect, Vect, Vect, Vect); 

    //function for reading text file
    char*  nextword(char*, char**);
    //delete atoms that are too close to each other
    void   DeleteAtom(int);
    void   DeleteAtomTe(int,int);
    void   DeleteAtomHeli(int);

    void   AddAtom(int);
    int    IsInsideLoop(int,Vect);
    void   ApplyPBC(Atom&);
    void   PBCShift(Vect&, Vect);
    double Fraction2Decimal(char*);
    int    Orthogonal(Vect, Vect, Vect);
    int    Right_handed(Vect, Vect, Vect);

    //check if it is a numerical parameter
    int    IsNumPar(char*);
    //find the nearest atom to the bottom face center of a SF tetrahedron 
    void   FindNearestAtom();

		void   ReadLatConst(int,char**,int,int);
		void   InitDirect(Direct&);

		//-------------------------------------------------------
		//to determine if atoms are shifted by a small distance
		//so that atoms are not perfectly align with dislocation
		//core
		//-------------------------------------------------------
		int atomShiftFlag;

		//-------------------------------------------------------
		//Miller or Miller-Baravis index
		//-------------------------------------------------------
		int    indexType; 
		void   ReadIndex(int,Direct&,char**,int,int);
		void   IndexTransform(Direct&);
		
    //-------------------------------------------------------
    //paramaters for dislocaiton loop
    //-------------------------------------------------------
    int    nloop; 
    Vect   *cdisl;       //center of dislocation loop
    double *rdisl;       //radius
    Direct *normdisl;    //normal
    Direct *lx;          //local x-axis on the slip plane
    Direct *ly;          //local y-axis on the slip plane 
    int    *ndisl;       //n segments for each loop
    Vect   **vdisl;      //vertex of discretized loop
    Direct *bv;          //Burgers vector
    double *bs;          //scaling factor of B-vector
    //paramaters for a rectangle dislocation loop 
    double *eratio;      // edge ratio dy/dx, default value is 1

    //-------------------------------------------------------
    //paramaters for SF tetrahedron
    //-------------------------------------------------------
    int     ntesf;
    Vect    *ctesf;
    double  *rtesf;
    Direct  *normtesf;
    Direct  *lxt;
    Direct  *lyt;
    Vect    **vtesf;

    //-------------------------------------------------------
    //paramaters for helical dislocation  
    //-------------------------------------------------------
    int      nheli;
    int      *ngridh;
    Triangle **gridh; 
    char     **fheli; //files that stores the mesh grid
    double   *bsh;
    //Vect   *bvh;
	Direct   *bvh;

    //-------------------------------------------------------
    //atoms in perfect crystals
    //-------------------------------------------------------
    Region regInfo;
    int    natom;
    int    ntype;
    double *mass;
    Atom   *atom;
    Atom   *atom_init;
    Vect   *u;

    //-------------------------------------------------------
    //input paramater labels: to very if all required parameters are provided
    //-------------------------------------------------------
    RegionParaLabel regParaLabel;
    LoopParaLabel   *loopParaLabel;
    HelixParaLabel  *helixParaLabel;
    TetraParaLabel  *tetraParaLabel;

    //size of the atomic region
    Direct x, y, z;           //direction of atomic region along x/y/z-axis
    double lat_a, lat_c;      //lattice constant
    int    pbcx, pbcy, pbcz;  //pbc flags: '1': pbc; '0': non-pbc
	int    resetpbcx, resetpbcy, resetpbcz;  //pbc is reset or not for output datafile
											 //flags: '1': reset to pbc, '0': no reset

	//matrial property
	double nu; // Poisson's ratio

    //parameter for deleting atoms that are too closely spaced 
    int    natom_del;
    double rCut;              //minimium seperation between atomic pairs
    int    *cellList;
    int    cellListLen;
    Nebr   *nebr;
    int    nebrListLen;

};

#endif





