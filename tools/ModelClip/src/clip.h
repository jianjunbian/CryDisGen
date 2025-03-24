/*----------------------------------------------

  Read atomic model, clip the model, then output
	the model;

	Data format: LAMMPS data file

	-----------------------------------------------*/

#ifndef CLIP_H
#define CLIP_H

#include "vect.h"

class Clip {
	public:
		Clip(int, char**);
		void ReadConfig(char *);
		void ClipModel();
		void Output(char *);

	private:
		char *nextword(char*, char**);

		//cylinder region
		char *axis;
		double coord1, coord2;
		double R;

		//origional model
		int ntype;
		double *mass;

		int natom;
		Atom *atom;
		Region regInfo;

		//clipped model 
		int natomClip;
		Atom *atomClip;
		Region regInfoClip;
};

#endif




