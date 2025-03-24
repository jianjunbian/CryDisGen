/*--------------------------------------------



	--------------------------------------------*/

#ifndef VECT_H
#define VECT_H

typedef struct {
	double x, y, z;
}Vect;

typedef struct {
	int id;
	int type;
	Vect coord;
	int stat;
}Atom;

typedef struct {
	double xlo, xhi, xprd;
	double ylo, yhi, yprd;
	double zlo, zhi, zprd;
}Region;

#endif
