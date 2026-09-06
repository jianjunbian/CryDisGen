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

#include "extra_math.h"
#include <math.h>

/*----------------------------------------------------*/
Vect VCross(Vect a, Vect b)
{
	Vect c; 
	c.x = a.y * b.z - a.z * b.y;
	c.y = a.z * b.x - a.x * b.z;
	c.z = a.x * b.y - a.y * b.x;
	return c;
}
/*----------------------------------------------------*/
double VDot(Vect a, Vect b)
{
	double s;
	s = a.x * b.x + a.y * b.y + a.z * b.z;
	return s;
}
/*----------------------------------------------------*/
double VLen(Vect v)
{
  double l;
	l = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	return l;
}
/*----------------------------------------------------*/
Vect VNorm(Vect v)
{
	double l = VLen(v);
	Vect vn;
	vn.x = v.x/l;
	vn.y = v.y/l;
	vn.z = v.z/l;
	return vn;
}
/*----------------------------------------------------*/
Vect VSub(Vect a, Vect b)
{
	Vect c; 
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	return c;
}
/*----------------------------------------------------*/
Vect VMul(Vect v, double s)
{
	Vect c;
	c.x = v.x * s;
	c.y = v.y * s;
	c.z = v.z * s;
	return c;
}
/*----------------------------------------------------*/
Vect VAdd(Vect a, Vect b)
{
	Vect c; 
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;
	return c;
}
/*----------------------------------------------------*/
Vect VMid(Vect a, Vect b)
{
	Vect c;
	c.x = (a.x + b.x)/2.;
	c.y = (a.y + b.y)/2.;
	c.z = (a.z + b.z)/2.;
	return c;
}
/*----------------------------------------------------*/
void NormalizeV(Vect &v)
{
	double l = VLen(v);
	v.x /= l;
	v.y /= l;
	v.z /= l;
}
/*----------------------------------------------------*/
void InitV(Vect &v)
{
	v.x = 0.;
	v.y = 0.;
	v.z = 0.;
}
/*----------------------------------------------------*/
Vect RotateX(Vect v, double alpha)
{
	double nv[3], ov[3], m[3][3];
	for(int i = 0; i < 3; i ++){
		for(int j = 0; j < 3; j ++){
			m[i][j] = 0.;
		}
	}
	for(int i = 0; i < 3; i ++){
		nv[i] = 0; 
		ov[i] = 0;
	}

	m[0][0] = 1;
	m[1][1] = cos(alpha); m[1][2] = -sin(alpha);
	m[2][1] = sin(alpha); m[2][2] = cos(alpha);

	ov[0] = v.x; ov[1] = v.y; ov[2] = v.z;
	for(int i = 0; i < 3; i ++){
		for(int j = 0; j < 3; j ++){
			nv[i] += m[i][j]*ov[j];
		}
	}
  Vect rv;
	rv.x = nv[0]; rv.y = nv[1]; rv.z = nv[2];
	return rv;
}
/*----------------------------------------------------*/
Vect RotateY(Vect v, double beta)
{
	double nv[3], ov[3], m[3][3];
	for(int i = 0; i < 3; i ++){
		for(int j = 0; j < 3; j ++){
			m[i][j] = 0;
		}
	}
	for(int i = 0; i < 3; i ++){
		ov[i] = 0; 
		nv[i] = 0;
	}
	m[0][0] = cos(beta); m[0][2] = sin(beta);
	m[1][1] = 1;
	m[2][0] = -sin(beta); m[2][2] = cos(beta);
	ov[0] = v.x; ov[1] = v.y; ov[2] = v.z;

	for(int i = 0; i < 3; i ++){
		for(int j = 0; j < 3; j ++){
			nv[i] += m[i][j]*ov[j];
		}
	}
	Vect rv;
	rv.x = nv[0]; rv.y = nv[1]; rv.z = nv[2];
	return rv;
}
/*----------------------------------------------------*/
Vect RotateZ(Vect v, double alpha)
{
	double nv[3], ov[3], m[3][3];
	for(int i = 0; i < 3; i ++){
		for(int j = 0; j < 3; j ++){
			m[i][j] = 0.;
		}
	}
	for(int i = 0; i < 3; i ++){
		nv[i] = 0; 
		ov[i] = 0;
	}

	m[0][0] = cos(alpha); m[0][1] = -sin(alpha);
	m[1][0] = sin(alpha); m[1][1] = cos(alpha);
	m[2][2] = 1;

	ov[0] = v.x; ov[1] = v.y; ov[2] = v.z;
	for(int i = 0; i < 3; i ++){
		for(int j = 0; j < 3; j ++){
			nv[i] += m[i][j]*ov[j];
		}
	}
  Vect rv;
	rv.x = nv[0]; rv.y = nv[1]; rv.z = nv[2];
	return rv;
}
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/
/*----------------------------------------------------*/

