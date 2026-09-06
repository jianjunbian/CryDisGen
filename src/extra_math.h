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

#ifndef EXTRA_MATH_H
#define EXTRA_MATH_H

#include "vect.h"

const double PI = 3.141592653589793238462643383279502884197169939375;

Vect   VCross(Vect, Vect);
double VDot(Vect, Vect);
double VLen(Vect);
Vect   VNorm(Vect);
void   NormalizeV(Vect&);
void   InitV(Vect&);
Vect   VSub(Vect, Vect);
Vect   VMul(Vect, double);
Vect   VAdd(Vect, Vect);
Vect   VMid(Vect, Vect);

Vect RotateX(Vect, double);
Vect RotateY(Vect, double);
Vect RotateZ(Vect, double);

#endif 
