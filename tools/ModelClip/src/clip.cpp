/*-------------------------------------------------------


  -------------------------------------------------------*/

#include "clip.h"
#include "vect.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*------------------------------------------------------*/
Clip::Clip(int argc, char **argv)
{
	if(argc != 6){
		printf("Usage:\n");
		printf("  Clip config x/y/z coord1 coord2 R\n");
		printf("\n");
		printf("Options:\n");
		printf("  config: atomic configuration\n");
		printf("   x/y/z: axis of cylinder region\n");
		printf("  coord1: position of cylidner axis in the other 2 dimensions\n");
		printf("  coord2: position of cylidner axis in the other 2 dimensions\n");
		printf("       R: radiu of cylinder region\n");
		exit(1);
	}

	axis = new char [8];
	strcpy(axis, argv[2]);
	coord1 = atof(argv[3]);
	coord2 = atof(argv[4]);
	R      = atof(argv[5]);

	ReadConfig(argv[1]);
	ClipModel();
	char outfile[64] = "clip.lmp";
	Output(outfile);
}
/*------------------------------------------------------*/
void Clip::ReadConfig(char *file)
{
  char *line, *copy;
  line = new char[1024];
  FILE *fp; 
  char *first, *next, *ptr;
  int nitem, narg, rowNum = 0;
  char **items = new char* [32*sizeof(char*)]; // assume that 20 items at most each line

	ntype = 0;

  fp = fopen(file, "r");
	if(fp == NULL){
		printf("Cannot find input file: %s\n", file);
		exit(1);
	}
  while(1) {  
    nitem = 0;
    if(fgets(line, 1024, fp) == NULL) break;
    first = nextword(line, &next);
    if(first == NULL) continue; // blank line  
    rowNum ++;  // add row number by one 
    ptr = next;
    items[0] = first;
    while(ptr) {
      nitem ++;
      items[nitem] = nextword(ptr, &next);
      if(!items[nitem]) break;
      ptr = next;
    }    
    narg = nitem;

    //if(!strcmp(items[nitem], '\0')) narg = nitem;  
    //else narg = nitem + 1;

		if(rowNum == 2){
			natom = atoi(items[0]);
			atom = new Atom[natom];
		}
		if(rowNum == 3){
			ntype = atoi(items[0]);
			mass = new double[ntype];
		}
		if(rowNum == 4){
			regInfo.xlo = atof(items[0]);
			regInfo.xhi = atof(items[1]);
			regInfo.xprd = regInfo.xhi - regInfo.xlo;
		}
		if(rowNum == 5){
			regInfo.ylo = atof(items[0]);
			regInfo.yhi = atof(items[1]);
			regInfo.yprd = regInfo.yhi - regInfo.ylo;
		}
		if(rowNum == 6){
			regInfo.zlo = atof(items[0]);
			regInfo.zhi = atof(items[1]);
			regInfo.zprd = regInfo.zhi - regInfo.zlo;
		}
		if(rowNum > 7 && rowNum <= 7 + ntype){
			mass[rowNum - 8] = atof(items[1]);
		}
		if(rowNum > 8 + ntype){
			int i = rowNum - 8 - ntype - 1;
		  atom[i].id      = atoi(items[0]);	
		  atom[i].type    = atoi(items[1]);	
			atom[i].coord.x = atof(items[2]);
			atom[i].coord.y = atof(items[3]);
			atom[i].coord.z = atof(items[4]);
			atom[i].stat    = 1;
		}
	}
  fclose(fp);
}
/*------------------------------------------------------*/
char* Clip::nextword(char *str, char **next)
{
  char *start, *stop;
  start = &str[strspn(str, " \t\n\v\f\r")];
  if(*start == '\0') return NULL; // space line
  stop = &start[strcspn(start," \t\n\v\f\r")];
  if(*stop == '\0') *next = NULL;
  else *next = stop+1;
  *stop = '\0';
  return start;
}
/*------------------------------------------------------*/
void Clip::ClipModel()
{
  natomClip = 0;
	atomClip  = NULL;

	Atom *buf;

	if(strcmp(axis, "x") == 0){
		for(int i = 0; i < natom; i ++){
			double dr1 = atom[i].coord.y - coord1;
			double dr2 = atom[i].coord.z - coord2;
			double drr = dr1 * dr1 + dr2 * dr2;
			if(drr < R*R){
				natomClip ++;
				buf = (Atom*)realloc(atomClip, natomClip*sizeof(Atom));
				if(buf == NULL){
					printf("Memory allocation error\n");
					exit(1);
				}
				atomClip = buf;
				atomClip[natomClip - 1] = atom[i];
			}
		}
	}
	else if(strcmp(axis, "y") == 0){
		for(int i = 0; i < natom; i ++){
			double dr1 = atom[i].coord.x - coord1;
			double dr2 = atom[i].coord.z - coord2;
			double drr = dr1 * dr1 + dr2 * dr2;
			if(drr < R*R){
				natomClip ++;
				buf = (Atom*)realloc(atomClip, natomClip*sizeof(Atom));
				if(buf == NULL){
					printf("Memory allocation error\n");
					exit(1);
				}
				atomClip = buf;
				atomClip[natomClip - 1] = atom[i];
			}
		}
	}
	else if(strcmp(axis, "z") == 0){
		for(int i = 0; i < natom; i ++){
			double dr1 = atom[i].coord.x - coord1;
			double dr2 = atom[i].coord.y - coord2;
			double drr = dr1 * dr1 + dr2 * dr2;
			if(drr < R*R){
				natomClip ++;
				buf = (Atom*)realloc(atomClip, natomClip*sizeof(Atom));
				if(buf == NULL){
					printf("Memory allocation error\n");
					exit(1);
				}
				atomClip = buf;
				atomClip[natomClip - 1] = atom[i];
			}
		}
	}

}
/*------------------------------------------------------*/
void Clip::Output(char *file)
{
	// output in lammps data-file format
  FILE *fp;
	char *fname = new char [64];
	strcpy(fname, file);
  fp = fopen(fname, "w");
  fprintf(fp, "LAMMPS description\n");
  fprintf(fp, "\n");
  fprintf(fp, "%d atoms\n", natomClip);
  fprintf(fp, "%d atom types\n", ntype);
  fprintf(fp, "\n");
  fprintf(fp, "%f %f xlo xhi\n", regInfo.xlo, regInfo.xhi); 
  fprintf(fp, "%f %f ylo yhi\n", regInfo.ylo, regInfo.yhi); 
  fprintf(fp, "%f %f zlo zhi\n", regInfo.zlo, regInfo.zhi); 
  fprintf(fp, "\n");

  //output data file body
  fprintf(fp, "Masses\n");
  fprintf(fp, "\n");
	for(int i = 0; i < ntype; i ++){
		fprintf(fp, "%d %f\n", i+1, mass[i]);
	}
  fprintf(fp, "\n");
  //atom coordination
  fprintf(fp, "Atoms\n");
  fprintf(fp, "\n");
	int id = 1;
  for(int i = 0; i < natomClip; i ++){
    fprintf(fp, "%d %d %f %f %f\n", id, atomClip[i].type, 
				atomClip[i].coord.x, atomClip[i].coord.y, atomClip[i].coord.z);
		id ++;
  }
  fclose(fp);

	printf("---------------Output file------------------\n");
	printf("  %s: in lammps data-file format\n", fname);
	printf("--------------------------------------------\n");

}
/*------------------------------------------------------*/
