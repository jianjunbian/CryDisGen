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

#include "createdisl.h"
#include "extra_math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/*------------------------------------------------------*/
CreateDisl::CreateDisl(int argc, char **argv)
{
  if(argc != 3){
    printf("Version 2.2 (May 2026)\n");
    printf("USAGE: \n");
    printf("\n");
    printf("    CryDisGen <parameter_file> <crystal_file>\n"); 
    printf("\n");
    printf("Description:\n");
    printf("\n");
    printf("    To create dislocations in fcc/bcc/hcp crystal structures\n");
    printf("\n");
    printf("Options:\n");
    printf("\n");
    printf("    parameter_file: A text file containing the dislocation setups\n");
	printf("                    Online documentation: www.crydisgen.com\n");
    printf("\n");
    printf("      crystal_file: A text file specifying the crystal structure, formatted\n");
    printf("                    to match the input requirements of the read_data command\n");
    printf("                    in LAMMPS\n");
    exit(1);
  }
  ReadParaFile(argv[1]);
  ReadConfig(argv[2]);

  natom_del = 0;
  //dislocation loop
  for(int i = 0; i < nloop; i ++){
    DiscretizeLoop(i);
    ApplyDisplacement(i);
    if(VDot(normdisl[i].v, VMul(bv[i].v, bs[i])) > 0 && fabs(VDot(normdisl[i].v, bv[i].v)) > 1e-6){
      printf("Deleting atoms ...\n");
      DeleteAtom(i);
    }
    else if(VDot(normdisl[i].v, VMul(bv[i].v, bs[i])) < 0 && fabs(VDot(normdisl[i].v, bv[i].v)) > 1e-6){
      printf("Adding atoms ...\n");
      AddAtom(i);
    }
  }
  //stacking fault tetrahedron 
  if(ntesf > 0) FindNearestAtom(); 
  for(int i = 0; i < ntesf; i ++){
    DiscretizeLoopTe(i);
    for(int faceid = 0; faceid < 4; faceid ++){
      ApplyDisplacementTe(i, faceid);
      if(faceid == 0) DeleteAtomTe(i, faceid); //operation is only needed for Frank loop
    }
  }
  //helical dislocation loop
  for(int i = 0; i < nheli; i ++){
    ReadGrid(i);
    ApplyDisplacementHeli(i);
    DeleteAtomHeli(i);
  }
  //output dislocation structure
  char *file = new char [64];
  strcpy(file, "outfile");
  Output(file);
}
/*------------------------------------------------------*/
void CreateDisl::ReadParaFile(char *file)
{
  char *line, *copy;
  int MAX_LENGTH = 1024; 
  line = new char[MAX_LENGTH];
  copy = new char[MAX_LENGTH];
  FILE *fp; 
  char *first, *next;  
  char *ptr;
  int nitem, narg, rowNum = 0;
  char **items = new char* [32*sizeof(char*)]; // assume that 20 items at most each line

  fp = fopen(file, "r");
  if(fp == NULL){
    printf("CanNOT find input file: %s\n", file);
    exit(1);
  }
  //to check if these labels appear only once
  int numLabelRegion = 0;
  int numLabelnLoop  = 0;
  int numLabelnTetra = 0;
  int numLabelnHelix = 0;
  //to store the number of input dislocaiton loops and other instances
  int loopLabel = 0; //dislocation loop label
  int id;
  int tesfLabel = 0; //stacking fault tetehedron
  int idt; 
  int heliLabel = 0; //helical dislocation
  int idh;
  //initialize input parameter label
  regParaLabel.region = 0;
  regParaLabel.x      = regParaLabel.y   = regParaLabel.y = 0;
  regParaLabel.lat    = regParaLabel.pbc = 0;
	//defautl: no resetting of PBC
	resetpbcx = resetpbcy = resetpbcz = 0;
  //default: there is no dislocation, stacking fault tetrahedron or helical dislocation
  nloop = 0;
  ntesf = 0;
  nheli = 0;
	//default value of Poisson's ratio
	nu = 0.; // to determine the elastic displacement is included (0) or not (!0) 
	//defautl flag for shifting atoms when calcualting elastic atomic displacement
	atomShiftFlag = 0;
  //-----Miller or Miller-Baravais index------------------------
  indexType = 0; //1: Miller index; 2: Miller-Baravais index
  //------------------------------------------------------------
  while(1) {  
    nitem = 0;
    if(fgets(line, MAX_LENGTH, fp) == NULL) { 
      break;
    } 
    //add row number by one 
    rowNum ++;                           
    //delete comment parts beginning with  "#", and store the rest in array 'copy'
    for(int i = 0; i < MAX_LENGTH; i ++) copy[i]='\0'; //initialize char array buffer
    int ncopy = strcspn(line, "#");
    strncpy(copy, line, ncopy);

    first = nextword(copy, &next);
    if(first == NULL) continue;           // blank line  
    ptr = next;
    items[0] = first;
    while(ptr) {
      nitem ++;
      items[nitem] = nextword(ptr, &next);
      if(!items[nitem]) break;
      ptr = next;
    }    
    narg = nitem;
    if(strcmp(items[0], "Region") == 0 || strcmp(items[0], "region") == 0) {
      if(nitem != 1){
        printf("Error: command 'Region' incorrect in line %d\n", rowNum);
        printf("       redundant parameter\n");
        exit(1);
      }
      numLabelRegion ++;
      if(numLabelRegion > 1){
        printf("Error: multiple 'Region' labels are encountered in line %d, only one is required\n", rowNum);
        exit(1);
      }
      regParaLabel.region = 1;
    }
    //read region info
    else if(strcmp(items[0], "-x") == 0){
      if( nitem != 4 && nitem != 5){
        printf("Error: -x parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      regParaLabel.x = 1;
      if(indexType == 0){                   //check Miller/Miller-Bravais indexes
        if(nitem == 4)      indexType = 1;  //Miller index
        else if(nitem == 5) indexType = 2;  //Miller-Baravais index
      }
      ReadIndex(indexType, x, items, nitem, rowNum);
    }
    else if(strcmp(items[0], "-y") == 0){
      if( nitem != 4 && nitem != 5){
        printf("Error: -y parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      regParaLabel.y = 1;
      if(indexType == 0){                   //check Miller/Miller-Bravais indexes
        if(nitem == 4)      indexType = 1;  //Miller index
        else if(nitem == 5) indexType = 2;  //Miller-Baravais index
      }
      ReadIndex(indexType, y, items, nitem, rowNum);
    }
    else if(strcmp(items[0], "-z") == 0){
      if(nitem != 4 && nitem != 5){
        printf("Error: -z parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      regParaLabel.z = 1;
      if(indexType == 0){                   //check Miller/Miller-Bravais indexes
        if(nitem == 4)      indexType = 1;  //Miller index
        else if(nitem == 5) indexType = 2;  //Miller-Baravais index
      }
      ReadIndex(indexType, z, items, nitem, rowNum);
    }
    else if(strcmp(items[0], "-lat") == 0){
      if(nitem != 2 && nitem != 3){
        printf("Error: -lat parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      regParaLabel.lat = 1;
      if(indexType == 0){
        if(nitem == 2)      indexType = 1;
        else if(nitem == 3) indexType = 2;
      }
      ReadLatConst(indexType, items, nitem, rowNum);
    }
    else if(strcmp(items[0], "-pbc") == 0){
      if(nitem != 4){
        printf("Error: -pbc parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      regParaLabel.pbc = 1; 
      if(!IsNumPar(items[1]) || !IsNumPar(items[2]) || !IsNumPar(items[3])){
        printf("non-numeric parameter appears in line %d\n", rowNum);
        exit(1);
      }
      pbcx = atoi(items[1]);
      pbcy = atoi(items[2]);
      pbcz = atoi(items[3]);
    }
		else if(strcmp(items[0], "-resetpbc") == 0){
			if(nitem != 4){
        printf("Error: -resetpbc parameter list incorrect in line %d\n", rowNum);
        exit(1);
			}
      if(!IsNumPar(items[1]) || !IsNumPar(items[2]) || !IsNumPar(items[3])){
        printf("non-numeric parameter appears in line %d\n", rowNum);
        exit(1);
      }
      resetpbcx = atoi(items[1]);
      resetpbcy = atoi(items[2]);
      resetpbcz = atoi(items[3]);
		}
		//read Poisson's ratio
		else if(strcmp(items[0], "-nu") == 0 || strcmp(items[0], "-Nu") == 0){
      if(nitem != 2 && nitem != 3){
        printf("Error: -nu parameter list incorrect in line %d\n", rowNum);
        exit(1);
			}
			if(!IsNumPar(items[1])){
				printf("non-numeric parameter encountered for Poisson's ratio in line %d\n", rowNum);
				exit(1);
			}
			nu = atof(items[1]); 
			if(nitem == 3){
				if(strcmp(items[2], "shift") != 0 && strcmp(items[2], "Shift") != 0){
					printf("Error: unkonw option \'%s\' in line %d, ", items[2], rowNum);
					printf("it should be 'shift' or 'Shift' here\n");
					exit(1);
				}
				atomShiftFlag = 1;
			}
		}
    //read the number of dislocation loop
    else if(strcmp(items[0], "nLoop") == 0 || strcmp(items[0], "nloop") == 0){
      if(nitem != 2){
        printf("Error: nLoop/nloop parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      numLabelnLoop ++;
      if(numLabelnLoop > 1){
        printf("Error: multiple 'nLoop' labels are encountered in line %d, only one is required\n",
            rowNum);
        exit(1);
      }
      nloop = atof(items[1]);
      if(nloop <= 0) {
        printf("Error: number of dislocation loop is NOT a positive integer in line %d\n", rowNum);
        exit(1);
      }
      //allocate memory
      cdisl    = new Vect   [nloop];
      rdisl    = new double [nloop];
      normdisl = new Direct [nloop];
      ndisl    = new int    [nloop];
      vdisl    = new Vect*  [nloop];
      lx       = new Direct [nloop];
      ly       = new Direct [nloop];
      bs       = new double [nloop];
      bv       = new Direct [nloop];
      eratio   = new double [nloop];
      //initialize Direct
      for(int i = 0; i < nloop; i ++){
        InitDirect(normdisl[i]);  
        InitDirect(lx[i]);  
        InitDirect(ly[i]);  
        InitDirect(bv[i]);  
      }
      //set default raito of two neighboring edges of a rectangular loop: ly/lx
			//so that it has a square shape
      for(int i = 0; i < nloop; i ++) eratio[i] = 1.;
      //labels for input parameters of dislocation loops     
      loopParaLabel = new LoopParaLabel [nloop];
      for(int i = 0; i < nloop; i ++){
        loopParaLabel[i].bv = loopParaLabel[i].bs = 0;
        loopParaLabel[i].dnorm = 0;
        loopParaLabel[i].dc = loopParaLabel[i].dr = 0;
        loopParaLabel[i].lx = loopParaLabel[i].ly = 0;
        loopParaLabel[i].ndisl = 0;
      }
    }
    else if(strcmp(items[0], "Loop") == 0 || strcmp(items[0], "loop") == 0) {
      if(nitem != 1){
        printf("Error: command 'Loop' incorrect in line %d\n", rowNum);
        printf("       redundant parameter\n");
        exit(1);
      }
      loopLabel ++;
      if(loopLabel > nloop){
        printf("----------------Error in parameter file-----------------\n");
        printf("  Number of dislocaiton loop set in 'nLoop N' \n");
        printf("  is different from that in following 'Loop' sections\n");
        printf("--------------------------------------------------------\n");
        exit(1);
      } 
      id = loopLabel - 1;
    }
    else if(strcmp(items[0], "-bv") == 0){
      if(nitem != 4 && nitem != 5){
        printf("Error: -bv parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      loopParaLabel[id].bv = 1;
      if(indexType == 0){                   //check Miller/Miller-Bravais indexes
        if(nitem == 4)      indexType = 1;  //Miller index
        else if(nitem == 5) indexType = 2;  //Miller-Baravais index
      }
      ReadIndex(indexType, bv[id], items, nitem, rowNum);
    }
    else if(strcmp(items[0], "-bs") == 0){
      if(nitem != 2){
        printf("Error: -bs parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      loopParaLabel[id].bs = 1;
      if(!IsNumPar(items[1])){
        printf("Error: non-numeric parameter appears in line %d\n", rowNum);
        exit(1);
      }
      bs[id] = Fraction2Decimal(items[1]);
      if(bs[id] == 0){
        printf("Error: zero value is NOT allowed in line %d\n", rowNum);
        exit(1);
      }
    }
    else if(strcmp(items[0], "-dnorm") == 0){
      if(nitem !=4 && nitem != 5){
        printf("Error: -dnorm parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      loopParaLabel[id].dnorm = 1;
      if(indexType == 0){                   //check Miller/Miller-Bravais indexes
        if(nitem == 4)      indexType = 1;  //Miller index
        else if(nitem == 5) indexType = 2;  //Miller-Baravais index
      }
      ReadIndex(indexType, normdisl[id], items, nitem, rowNum);
    }
    else if(strcmp(items[0], "-dc") == 0){
      if(nitem != 4){
        printf("Error: -dc parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      loopParaLabel[id].dc = 1;
      if(!IsNumPar(items[1]) || !IsNumPar(items[2]) || !IsNumPar(items[3])){
        printf("Error: non-numeric parameter appears in line %d\n", rowNum);
        exit(1);
      }
      cdisl[id].x = atof(items[1]);
      cdisl[id].y = atof(items[2]);
      cdisl[id].z = atof(items[3]);
    }
    else if(strcmp(items[0], "-dr") == 0){
      if(nitem != 2){
        printf("Error: -dr parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      loopParaLabel[id].dr = 1;
      if(!IsNumPar(items[1])){
        printf("Error: non-numeric parameter appears in line %d\n", rowNum);
        exit(1);
      }
      rdisl[id] = atof(items[1]);
    }
    else if(strcmp(items[0], "-lx") == 0){
      if(nitem != 4 && nitem != 5 && nitem != 2){
        printf("Error: -lx parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      loopParaLabel[id].lx = 1;
      if(nitem != 2) { 
        if(indexType == 0){                   //check Miller/Miller-Bravais indexes
          if(nitem == 4)      indexType = 1;  //Miller index
          else if(nitem == 5) indexType = 2;  //Miller-Baravais index
        }
        ReadIndex(indexType, lx[id], items, nitem, rowNum);
      } 
      else if(nitem == 2 ){
        //vector '-lx' will be calculated by CryDisGen in DiscretizeLoop(int)
        if(strcmp(items[1], "null") != 0 && strcmp(items[1], "Null") != 0){
          printf("Error: -lx parameter incorrect in line %d\n", rowNum);
          exit(1);
        }
      }
    }
    else if(strcmp(items[0], "-ly") == 0){
      if(nitem != 4 && nitem != 5 && nitem != 2){
        printf("Error: -ly parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      loopParaLabel[id].ly = 1;
      if(nitem != 2){
        if(indexType == 0){                   //check Miller/Miller-Bravais indexes
          if(nitem == 4)      indexType = 1;  //Miller index
          else if(nitem == 5) indexType = 2;  //Miller-Baravais index
        }
        ReadIndex(indexType, ly[id], items, nitem, rowNum);
      }
      else if(nitem == 2 ){
        //vector '-ly' will be calculated by CryDisGen
        if(strcmp(items[1], "null") != 0 && strcmp(items[1], "Null") != 0){
          printf("Error: -ly parameter incorrect in line %d\n", rowNum);
          exit(1);
        }
      }
    }
    else if(strcmp(items[0], "-ndisl") == 0){
      if(nitem != 2 && nitem != 3){
        printf("Error: -ndisl parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      loopParaLabel[id].ndisl = 1;
      if(!IsNumPar(items[1])){
        printf("Error: non-numeric parameter appears in line %d\n", rowNum);
        exit(1);
      }
      ndisl[id] = atoi(items[1]);
      vdisl[id] = new Vect[ndisl[id]];
      //for sqare-shaped loop, the second parameter is used to scale the edge, dy/dx 
      if(nitem == 3 ){
        if(!IsNumPar(items[2])){
          printf("Error: non-numeric parameter appears in line %d\n", rowNum);
          exit(1);
        }
        if(ndisl[id] == 4) //eratio is only used for rectangular loop
          eratio[id] = atof(items[2]);
        else {
          printf("Error: -ndisl: the ratio parameter is only used for square-shaped loop (# of segment = 4)\n ");
          exit(1);
        }
      }
    }
    // stacking fault tetrahedron
    else if(strcmp(items[0], "nTetra") == 0 || strcmp(items[0], "ntetra")==0 || 
        strcmp(items[0], "nTetrahedron") == 0 || strcmp(items[0], "ntetrahedron") == 0){
      if(indexType != 1) {
        printf("HCP lattice doesn't support 'nTetrahedron' tag\n");
        exit(1);
      }
      if(nitem != 2){
        printf("Error: 'nTetra/ntetra' parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      numLabelnTetra ++;
      if(numLabelnTetra > 1){
        printf("Error: multiple 'nTetrahedron/nTetra' labels are encountered in line %d, only one is required\n", rowNum);
        exit(1);
      }
      ntesf = atoi(items[1]);
      if(ntesf <= 0){
        printf("----------------------------Error in parameter file----------------------------\n");
        printf("Error: the number of stacking fault tetrahedron is NOT a positive interger.\n");
        printf("-------------------------------------------------------------------------------\n");
        exit(1);
      }
      ctesf    = new Vect   [ntesf];
      rtesf    = new double [ntesf];
      normtesf = new Direct [ntesf];
      lxt      = new Direct [ntesf];
      lyt      = new Direct [ntesf];
      vtesf    = new Vect*  [ntesf];  // vertices
      for(int i = 0; i < ntesf; i ++){
        InitDirect(normtesf[i]);
        InitDirect(lxt[i]);
        InitDirect(lyt[i]);
      }
      for(int i = 0; i < ntesf; i ++){
        vtesf[i] = new Vect [4];
      }
      //input parameter label
      tetraParaLabel = new TetraParaLabel[ntesf];
      for(int i = 0; i < ntesf; i ++){
        tetraParaLabel[i].tc = tetraParaLabel[i].tr = 0;
        tetraParaLabel[i].tnorm = 0;
        tetraParaLabel[i].lxt = tetraParaLabel[i].lyt = 0;
      }
    }
    else if(strcmp(items[0], "Tetra") == 0 || strcmp(items[0], "tetra") == 0 ||
        strcmp(items[0], "Tetrahedron") == 0 || strcmp(items[0], "tetrahedron") == 0) {
      if(indexType != 1) {
        printf("HCP lattice doesn't support '%s' tag\n", items[0]);
        //printf("HCP lattice doesn't support 'nTetrahedron' tag\n");
        exit(1);
      }
      if(nitem != 1){
        printf("Error: command 'Tetrahedron/Tetra' incorrect in line %d\n", rowNum);
        printf("       redundant parameter\n");
        exit(1);
      }
      tesfLabel ++;
      if(tesfLabel > ntesf){
        printf("----------------Error in parameter file-----------------\n");
        printf("  Number of stacking fault tetrahedron %d in 'nTetra N' \n", ntesf);
        printf("  is different from %d that in following 'Tetra' sections\n",tesfLabel);
        printf("--------------------------------------------------------\n");
        exit(1);
      } 
      idt = tesfLabel - 1;
    }
    else if(!strcmp(items[0], "-tc")){
      if(indexType != 1) {
        printf("HCP lattice doesn't support '%s' tag\n", items[0]);
        exit(1);
      }
      if(nitem != 4){
        printf("Error: -tc parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      tetraParaLabel[idt].tc = 1;
      if(!IsNumPar(items[1]) || !IsNumPar(items[2]) || !IsNumPar(items[3])){
        printf("Error: non-numeric parameter appears in line %d\n", rowNum);
        exit(1);
      }
      ctesf[idt].x = atof(items[1]);
      ctesf[idt].y = atof(items[2]);
      ctesf[idt].z = atof(items[3]);
    }
    else if(!strcmp(items[0], "-tr")){
      if(indexType != 1) {
        printf("HCP lattice doesn't support '%s' tag\n", items[0]);
        //printf("HCP lattice doesn't support 'nTetrahedron' tag\n");
        exit(1);
      }
      if(nitem != 2){
        printf("Error: -tr parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      tetraParaLabel[idt].tr = 1;
      if(!IsNumPar(items[1])){
        printf("Error: non-numeric parameter appears in line %d\n", rowNum);
        exit(1);
      }
      rtesf[idt] = atof(items[1]);
			if(rtesf[idt] <= 0){
        printf("Error: invalid size of stackign fault tetrahedron in line %d\n", rowNum);
				exit(1);
			}
    }
    else if(!strcmp(items[0], "-tnorm")){
      if(indexType != 1) {
        printf("HCP lattice doesn't support '%s' tag\n", items[0]);
        exit(1);
      }
      if(nitem != 4){
        printf("Error: -tnorm parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      tetraParaLabel[idt].tnorm = 1;

      if(indexType == 0){                   //check Miller/Miller-Bravais indexes
        if(nitem == 4)      indexType = 1;  //Miller index
        else if(nitem == 5) indexType = 2;  //Miller-Baravais index
      }
      ReadIndex(indexType, normtesf[idt], items, nitem, rowNum);
    }
    else if(!strcmp(items[0], "-lxt")){
      if(indexType != 1) {
        printf("HCP lattice doesn't support '%s' tag\n", items[0]);
        exit(1);
      }
      if(nitem != 4){
        printf("Error: -lxt parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      tetraParaLabel[idt].lxt = 1;
      if(indexType == 0){                   //check Miller/Miller-Bravais indexes
        if(nitem == 4)      indexType = 1;  //Miller index
        else if(nitem == 5) indexType = 2;  //Miller-Baravais index
      }
      ReadIndex(indexType, lxt[idt], items, nitem, rowNum);
    }
    else if(!strcmp(items[0], "-lyt")){
      if(indexType != 1) {
        printf("HCP lattice doesn't support '%s' tag\n", items[0]);
        exit(1);
      }
      if(nitem != 4){
        printf("Error: -lyt parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      tetraParaLabel[idt].lyt = 1;
      if(indexType == 0){                   //check Miller/Miller-Bravais indexes
        if(nitem == 4)      indexType = 1;  //Miller index
        else if(nitem == 5) indexType = 2;  //Miller-Baravais index
      }
      ReadIndex(indexType, lyt[idt], items, nitem, rowNum);
    }
    else if(!strcmp(items[0], "nHelix") || !strcmp(items[0], "nhelix")){
      if(nitem != 2){
        printf("Error: nHelix/nhelix parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      numLabelnHelix ++;
      if(numLabelnHelix > 1){
        printf("Error: multiple 'nHelix' labels are encountered in line %d, only one is required\n", 
            rowNum);
        exit(1);
      }
      nheli   = atoi(items[1]);
      if(nheli <= 0){
        printf("----------------Error in parameter file-----------------\n");
        printf("Error: number of dislocation helix is NOT a positive number in line %d\n", rowNum);
        printf("--------------------------------------------------------\n");
        exit(1);
      }
      ngridh  = new int [nheli];
      gridh   = new Triangle *[nheli];
      fheli   = new char*[nheli];
      for(int i = 0; i < nheli; i ++){ fheli[i] = new char [64]; }
      bsh     = new double [nheli];
      bvh     = new Direct [nheli];
      //input parameter label
      helixParaLabel = new HelixParaLabel[nheli];
      for(int i = 0; i < nheli; i ++){
        helixParaLabel[i].bvh = helixParaLabel[i].bsh = 0;
        helixParaLabel[i].f = 0;
      }
    }
    else if(!strcmp(items[0], "Helix") || !strcmp(items[0], "helix")) {
      if(nitem != 1){
        printf("Error: command 'Helix/helix' incorrect in line %d\n", rowNum);
        printf("       redundant parameter\n");
        exit(1);
      }
      heliLabel ++;
      if(heliLabel > nheli){
        printf("----------------Error in parameter file-----------------\n");
        printf("  Number of helical dislocaiton set in cmd 'nHelical N' \n");
        printf("  is different from that in following 'Helical' sections\n");
        printf("--------------------------------------------------------\n");
        exit(1);
      } 
      idh = heliLabel - 1;
    }
    else if(!strcmp(items[0], "-bvh")){
      if(nitem != 4 && nitem != 5){
        printf("Error: -bvh parameter list incorrect\n");
        exit(1);
      }
      helixParaLabel[idh].bvh = 1;
      if(indexType == 0){                   //check Miller/Miller-Bravais indexes
        if(nitem == 4)      indexType = 1;  //Miller index
        else if(nitem == 5) indexType = 2;  //Miller-Baravais index
      }
      ReadIndex(indexType, bvh[idh], items, nitem, rowNum);
    }
    else if(!strcmp(items[0], "-bsh")){
      if(nitem != 2){
        printf("Error: -bsh parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      helixParaLabel[idh].bsh = 1;
      if(!IsNumPar(items[1])){
        printf("Error: non-numeric parameter appears in line %d\n", rowNum);
        exit(1);
      }
      bsh[idh] = Fraction2Decimal(items[1]);
      if(bsh[id] == 0){
        printf("Error: zero value is NOT allowed in line %d\n", rowNum); exit(1);
      }
    }
    else if(!strcmp(items[0], "-f")){
      if(nitem != 2){
        printf("Error: -f parameter list incorrect in line %d\n", rowNum);
        exit(1);
      }
      helixParaLabel[idh].f = 1;
      strcpy(fheli[idh], items[1]);
    }
    else {
      printf("Error: Unknown command \"%s\" in line %d in parameter file\n", items[0], rowNum );
      exit(1);
    }
  }
  fclose(fp); //finish reading the parameter file

  //to check the number of dislocation structure 
	if(nloop != loopLabel){
		printf("----------------Error in parameter file-----------------\n");
		printf("  'nLoop N' is NOT consistent with 'Loop' labels\n");
		printf("  'N' must equal to the total number of 'Loop' labels\n");
		printf("--------------------------------------------------------\n");
		exit(1);
	}
	if(ntesf != tesfLabel){
		printf("----------------Error in parameter file-----------------\n");
		printf("  'nTetrahedron N' is NOT consistent with 'Tetrahedron' labels\n");
		printf("  'N' must equal to the total number of 'Tetrahedron' labels\n");
		printf("--------------------------------------------------------\n");
		exit(1);
	}
	if(nheli != heliLabel){
		printf("----------------Error in parameter file-----------------\n");
		printf("  'nHelix N' is NOT consistent with 'Helix' labels\n");
		printf("  'N' must equal to the total number of 'Helix' labels\n");
		printf("--------------------------------------------------------\n");
		exit(1);
	}
	/*
  //to check the number of dislocation structure 
  if(loopLabel != nloop) {
    printf("Error: Number of tag \"Loop\" must equal to that set by \"nLoop n\"\n");
    exit(1);
  }
  if(tesfLabel != ntesf) {
    printf("Error: Number of tag \"Tetrahedron \" must equal to that set by \"nTetrahedron n\"\n");
    exit(1);
  }
  if(heliLabel != nheli) {
    printf("Error: Number of tag \"Helix\" must equal to that set by \"nHelix n\"\n");
    exit(1);
  }
	*/

  //transform Miller/Miller-Baravis index to vectors in Cartesian coordinate system
  IndexTransform(x);
  IndexTransform(y);
  IndexTransform(z);
  for(int i = 0; i < nloop; i ++){
    IndexTransform(normdisl[i]);
    IndexTransform(lx[i]);
    IndexTransform(ly[i]);
    IndexTransform(bv[i]);
  }
  for(int i = 0; i < ntesf; i ++){
    IndexTransform(lxt[i]);
    IndexTransform(lyt[i]);
    IndexTransform(normtesf[i]);
  }
  for(int i = 0; i < nheli; i ++){
    IndexTransform(bvh[i]);
  }

  //check if parameters of atomic region are read in
  if(regParaLabel.region != 1){
    printf("Atomic region is NOT set, 'Region' label is NOT found\n");
    exit(1);
  }
  else if(regParaLabel.x != 1){
    printf("x-axis of atomic region is NOT set\n");
    exit(1);
  }
  else if(regParaLabel.y != 1){
    printf("y-axis of atomic region is NOT set\n");
    exit(1);
  }
  else if(regParaLabel.z != 1){
    printf("z-axis of atomic region is NOT set\n");
    exit(1);
  }
	else if(regParaLabel.pbc != 1){
		printf("PBC of atomic region is NOT set\n");
		exit(1);
	}
  //check if the parameters of dislocation loops are all read in
  for(int i = 0; i < nloop; i ++){
    if(loopParaLabel[i].bv != 1 || loopParaLabel[i].bs != 1 || loopParaLabel[i].dnorm != 1 ||
        loopParaLabel[i].dc != 1 || loopParaLabel[i].dr != 1 ||
        loopParaLabel[i].lx != 1 || loopParaLabel[i].ly != 1 || loopParaLabel[i].ndisl != 1){
      printf("-------------------Error in parameter_file-----------------------------\n");
      printf("  Loop ID: %d\n", i + 1);
      printf("  The parameters of dislocation loop are NOT fully set\n");
      printf("  Full parameter list (1: set; 0: NOT set):\n");
      printf("    -bv: %d, -bs: %d, -dnorm: %d, -dc: %d, -dr: %d, -lx: %d, -ly: %d, -ndisl: %d\n ",
          loopParaLabel[i].bv, loopParaLabel[i].bs, loopParaLabel[i].dnorm, loopParaLabel[i].dc,
          loopParaLabel[i].dr, loopParaLabel[i].lx, loopParaLabel[i].ly,    loopParaLabel[i].ndisl);
      printf("-----------------------------------------------------------------------\n");
      exit(1);
    }
  }
  //check if the sparameters of stacking fault tetrahedrons are all read in
  for(int i = 0; i < ntesf; i ++){
    if(tetraParaLabel[i].tc != 1 || tetraParaLabel[i].tr != 1 || tetraParaLabel[i].tnorm != 1 ||
        tetraParaLabel[i].lxt != 1 || tetraParaLabel[i].lyt != 1 ){
      printf("-------------------Error in parameter_file-----------------------------\n");
      printf("  Stacking fault tetrahedron ID: %d\n", i + 1);
      printf("  The parameters of stacking fault tetrahedron are NOT fully set\n");
      printf("  Full parameter list (1: set; 0: NOT set):\n");
      printf("    -tc: %d, -tr: %d, -tnorm: %d, -lxt: %d, -lyt: %d\n", 
          tetraParaLabel[i].tc, tetraParaLabel[i].tr, tetraParaLabel[i].tnorm,
          tetraParaLabel[i].lxt, tetraParaLabel[i].lyt);
      printf("-----------------------------------------------------------------------\n");
      exit(1);
    }
  }
  //check if the parameters of helixes are all read in
  for(int i = 0; i < nheli; i ++){
    if(helixParaLabel[i].bvh != 1 || helixParaLabel[i].bsh != 1 || helixParaLabel[i].f != 1){ 
      printf("-------------------Error in parameter_file-----------------------------\n");
      printf("  Helix ID: %d\n", i + 1);
      printf("  The parameters of dislocation helixes are NOT fully set\n");
      printf("  Full parameter list (1: set; 0: NOT set):\n");
      printf("    -bvh: %d, -bsh: %d, -f: %d\n", 
          helixParaLabel[i].bvh, helixParaLabel[i].bsh, helixParaLabel[i].f);
      printf("-----------------------------------------------------------------------\n");
      exit(1);
    }
  }
  //check if the three axises of the crystal region are orthogonal
  if(!Orthogonal(x.v, y.v, z.v)){
      printf("-------------------Error in parameter_file----------------\n");
      printf("  Three axises of the crystal region are NOT orthogonal\n");
      printf("----------------------------------------------------------\n");
      exit(1);
  }
  //check if '-lx' and '-ly' both are set to null or if the axises of local 
  //coordinate system are muturally orthogonal
  Vect xn, yn, zn;
  for(int i = 0; i < nloop; i ++){
    if(VLen(lx[i].v)*VLen(ly[i].v) == 0){    //at least one equals to 0
      if(VLen(lx[i].v) + VLen(ly[i].v) == 0) //check if both equal to 0
        continue; 
      else{ 
        //if only one of '-lx' and '-ly' is set to null, print error info 
        printf("-------------------Error in parameter_file----------------\n");
        printf("  Dislocation loop ID: %d\n", i + 1);
        printf("  '-lx' and '-ly' must be set to 'Null/null' in pairs\n");
        printf("----------------------------------------------------------\n");
        exit(1);
      }
    }
    else {
      xn = lx[i].v;
      yn = ly[i].v;
      zn = normdisl[i].v;
      if(!Orthogonal(xn, yn, zn)){
        printf("-------------------Error in parameter_file----------------\n");
        printf("  Dislocation loop ID: %d\n", i + 1);
        printf("  Local coordinate axises of dislocaiton loop\n");
        printf("  are NOT orthogonal, i.e., dnorm, lx and ly\n");
        printf("----------------------------------------------------------\n");
        exit(1);
      }
    }
  }
  for(int i = 0; i < ntesf; i ++){
    xn = lxt[i].v;
    yn = lyt[i].v;
    zn = normtesf[i].v;
    if(!Orthogonal(xn, yn, zn)){
      printf("-------------------Error in parameter_file----------------\n");
      printf("  Stacking Fault Tetrahedron ID: %d\n", i + 1);
      printf("  Local coordinate axises of stacking fault tetrahedron\n");
      printf("  are NOT orthogonal, i.e., tnorm, lxt and lyt\n");
      printf("----------------------------------------------------------\n");
      exit(1);
    }
  }
  //check if the axises of the golbal and local coordinate systems are RIGHT handed
  if(!Right_handed(x.v, y.v, z.v)){
      printf("-------------------Error in parameter_file----------------\n");
      printf("  Three axises of the crystal region are NOT right handed\n");
      printf("----------------------------------------------------------\n");
      exit(1);
  }
  for(int i = 0; i < nloop; i ++){
    if(VLen(lx[i].v)*VLen(ly[i].v) == 0) continue;
    xn = lx[i].v;
    yn = ly[i].v;
    zn = normdisl[i].v;
    if(!Right_handed(xn, yn, zn)){
      printf("-------------------Error in parameter_file----------------\n");
      printf("  Dislocation loop ID: %d\n", i + 1);
      printf("  Local coordinate axises of dislocaiton loop\n");
      printf("  are NOT right handed, i.e., dnorm, lx and ly\n");
      printf("----------------------------------------------------------\n");
      exit(1);
    }
  }
  for(int i = 0; i < ntesf; i ++){
    xn = lxt[i].v;
    yn = lyt[i].v;
    zn = normtesf[i].v;
    if(!Right_handed(xn, yn, zn)){
      printf("-------------------Error in parameter_file----------------\n");
      printf("  Stacking Fault Tetrahedron ID: %d\n", i + 1);
      printf("  Local coordinate axises of stacking fault tetrahedron\n");
      printf("  are NOT right handed, i.e., tnorm, lxt and lyt\n");
      printf("----------------------------------------------------------\n");
      exit(1);
    }
  }
}
/*------------------------------------------------------*/
void CreateDisl::DiscretizeLoop(int loopid)
{
  int nv = ndisl[loopid];       
  double dtheta = 2*PI / nv;

  if(nv == 3){
    for(int i = 0; i < nv; i ++){
      vdisl[loopid][i].x = rdisl[loopid] * cos(dtheta*i - PI/6.);
      vdisl[loopid][i].y = rdisl[loopid] * sin(dtheta*i - PI/6.);
      vdisl[loopid][i].z = 0;
    }   
  }else if(nv == 4){
    for(int i = 0; i < nv; i ++){
      /*
      //square shape
      vdisl[loopid][i].x = rdisl[loopid] * cos(dtheta*i - PI/4.);
      vdisl[loopid][i].y = rdisl[loopid] * sin(dtheta*i - PI/4.);
      */
      vdisl[loopid][i].z = 0;
    }   
    //rectanngle when scaling factor is changed, dr is half size of edge length
    vdisl[loopid][0].x = rdisl[loopid]; vdisl[loopid][0].y =-eratio[loopid] * rdisl[loopid];
    vdisl[loopid][1].x = rdisl[loopid]; vdisl[loopid][1].y = eratio[loopid] * rdisl[loopid];
    vdisl[loopid][2].x =-rdisl[loopid]; vdisl[loopid][2].y = eratio[loopid] * rdisl[loopid];
    vdisl[loopid][3].x =-rdisl[loopid]; vdisl[loopid][3].y =-eratio[loopid] * rdisl[loopid];
  }
  else {
    for(int i = 0; i < nv; i ++){
      vdisl[loopid][i].x = rdisl[loopid] * cos(dtheta*i);
      vdisl[loopid][i].y = rdisl[loopid] * sin(dtheta*i);
      vdisl[loopid][i].z = 0;
    }   
  }

	//change the coordinate to global Coordinate system 
  Vect xn, yn, zn; //new axises 
  if(VLen(lx[loopid].v) != 0 && VLen(ly[loopid].v) != 0){
    zn = normdisl[loopid].v;
    yn =       ly[loopid].v;
    xn =       lx[loopid].v;
  }
  else{
    zn = normdisl[loopid].v;
    yn = VCross(zn, z.v);
    if(VLen(yn) < 1e-10) //in case yn is parallel to z.v
      yn = VCross(zn, x.v);
    xn = VCross(yn, zn);
  }

  double m[3][3], np[3], op[3];
  m[0][0] = VDot(xn, x.v)/VLen(xn)/VLen(x.v);
  m[1][0] = VDot(xn, y.v)/VLen(xn)/VLen(y.v);
  m[2][0] = VDot(xn, z.v)/VLen(xn)/VLen(z.v);

  m[0][1] = VDot(yn, x.v)/VLen(yn)/VLen(x.v);
  m[1][1] = VDot(yn, y.v)/VLen(yn)/VLen(y.v);
  m[2][1] = VDot(yn, z.v)/VLen(yn)/VLen(z.v);

  m[0][2] = VDot(zn, x.v)/VLen(zn)/VLen(x.v);
  m[1][2] = VDot(zn, y.v)/VLen(zn)/VLen(y.v);
  m[2][2] = VDot(zn, z.v)/VLen(zn)/VLen(z.v);

  for(int n = 0; n < ndisl[loopid]; n ++){
    for(int i = 0; i < 3; i ++) np[i] = 0; 
    op[0] = vdisl[loopid][n].x; 
    op[1] = vdisl[loopid][n].y; 
    op[2] = vdisl[loopid][n].z;
    for(int i = 0; i < 3; i ++){
      for(int j = 0; j < 3; j ++){
        np[i] += m[i][j] * op[j];
      }
    }
    vdisl[loopid][n].x = np[0]; 
    vdisl[loopid][n].y = np[1]; 
    vdisl[loopid][n].z = np[2];
    vdisl[loopid][n] = VAdd(vdisl[loopid][n], cdisl[loopid]);
  }
}
/*-------------------------------------------------------------*/
double CreateDisl::SolidAngleTetrahedron(Tetrahedron te)
{
  //top vertex v0;  bottom v1-v2-v3
  double alpha, beta, gama;
  Vect v10, v12, v13;
  Vect n102, n103;
  v10   = VSub(te.v0, te.v1);
  v12   = VSub(te.v2, te.v1);//
  v13   = VSub(te.v3, te.v1);
  n102  = VCross(v10, v12);
  n103  = VCross(v10, v13);
  alpha = acos( VDot(n102, n103) / VLen(n102) / VLen(n103) );

  Vect v20, v21, v23;
  Vect n201, n203;
  v20  = VSub(te.v0, te.v2);
  v21  = VSub(te.v1, te.v2);
  v23  = VSub(te.v3, te.v2);
  n201 = VCross(v20, v21);
  n203 = VCross(v20, v23);
  beta = acos( VDot(n201, n203) / VLen(n201) / VLen(n203) );

  Vect v30, v31, v32;
  Vect n301, n302;
  v30  = VSub(te.v0, te.v3);
  v31  = VSub(te.v1, te.v3);
  v32  = VSub(te.v2, te.v3);
  n301 = VCross(v30, v31);
  n302 = VCross(v30, v32);
  gama = acos( VDot(n301, n302) / VLen(n301) / VLen(n302) );
  
  double angle = alpha + beta + gama - PI;
  return angle;
}
/*------------------------------------------------------*/
double CreateDisl::SolidAngleTriangle(Triangle t, Vect p)
{
  double val, angle;
  double dist;
  Vect v0p, v1p, v2p;
  Vect v01, v10, v02, v20, v12, v21;
  double angle0, angle1, angle2;

  v0p = VSub(p, t.v0);
  v1p = VSub(p, t.v1);
  v2p = VSub(p, t.v2);

  v01 = VSub(t.v1, t.v0);
  v10 = VSub(t.v0, t.v1);
  v02 = VSub(t.v2, t.v0);
  v20 = VSub(t.v0, t.v2);
  v12 = VSub(t.v2, t.v1);
  v21 = VSub(t.v1, t.v2);

  angle0 = acos( VDot(v01, v02) / VLen(v01) / VLen(v02) );
  angle1 = acos( VDot(v10, v12) / VLen(v10) / VLen(v12) );
  angle2 = acos( VDot(v21, v20) / VLen(v21) / VLen(v20) );

  //case 1: p overlaps with one vertex
  if(VLen(v0p) < 1e-6){
    return angle0;
  }else if(VLen(v1p) < 1e-6){
    return angle1;
  }else if(VLen(v2p) < 1e-6){
    return angle2;
  }
  //case 2: p is at one of the edges
  if( (VDot(v01, v0p) > 0 && VDot(v10, v1p) > 0 && VLen(VCross(v01, v0p)) < 1e-6)  || 
      (VDot(v02, v0p) > 0 && VDot(v20, v2p) > 0 && VLen(VCross(v02, v0p)) < 1e-6)  ||
      (VDot(v12, v1p) > 0 && VDot(v21, v2p) > 0 && VLen(VCross(v12, v1p)) < 1e-6) )  {
    return PI;
  }
  //case 3: p is inside the triangle
  double angle10p, anglep02, anglep12, anglep10;
  angle10p = acos( VDot(v01, v0p) / VLen(v01) / VLen(v0p) );  
  anglep02 = acos( VDot(v02, v0p) / VLen(v02) / VLen(v0p) );
  anglep12 = acos( VDot(v12, v1p) / VLen(v12) / VLen(v1p) );
  anglep10 = acos( VDot(v1p, v10) / VLen(v1p) / VLen(v10) );
  if ( angle10p < angle0 && anglep02 < angle0 && 
      anglep10 < angle1 && anglep12 < angle1){
    return 2*PI;
  }
  return 0;
}
/*------------------------------------------------------*/
double CreateDisl::AtomSolidAngle(int loopid, Vect coord)
{
  //check if the atom is on plane of dislocaiton loop  
  Vect vdisl01 = VSub(vdisl[loopid][1], vdisl[loopid][0]);
  Vect vdisl02 = VSub(vdisl[loopid][2], vdisl[loopid][0]);
  Vect cross   = VCross(vdisl01, vdisl02);
  Vect vdisl0a = VSub(coord, vdisl[loopid][0]);
  //if the vector 'cross' is not normalized, it may produce a vector with 
  //a large magnitude and affect the value of variable 'val' 
  cross   = VNorm(cross);
  vdisl0a = VNorm(vdisl0a);
  double val = VDot(cross, vdisl0a);

  Tetrahedron te;
  Triangle t;
  double sAngle = 0;
  int next; 

  //case 1: atom on the dislocation plane
  if(fabs(val) < 1e-6){ 
    for(int i = 0; i < ndisl[loopid]; i ++){
      t.v0 = cdisl[loopid];
      t.v1 = vdisl[loopid][i];
      next = i + 1; 
      if(next == ndisl[loopid]) next = 0;
      t.v2 = vdisl[loopid][next];
      sAngle += SolidAngleTriangle(t, coord); 
    }
    return sAngle;
  }else{
  //case 2: atom lays out of the dislocation plane
    for(int i = 0; i < ndisl[loopid]; i ++){
      te.v0 = coord; 
      te.v1 = cdisl[loopid];
      te.v2 = vdisl[loopid][i];
      next = i + 1;
      if(next == ndisl[loopid]) next = 0;
      te.v3 = vdisl[loopid][next];
      sAngle += SolidAngleTetrahedron(te);
    }
    if(val > 0) return sAngle;
    else return -sAngle;
  }
}
/*------------------------------------------------------*/
void CreateDisl::ApplyDisplacement(int loopid)
{
  Vect b = VMul(bv[loopid].v, bs[loopid]);
  //transform Burgers vector to vector in Cartesian coordinate system
  double vlen = VLen(b);
  double cosbx = VDot(b, x.v)/vlen/VLen(x.v);
  double cosby = VDot(b, y.v)/vlen/VLen(y.v);
  double cosbz = VDot(b, z.v)/vlen/VLen(z.v);
  b.x = vlen * cosbx;
  b.y = vlen * cosby;
  b.z = vlen * cosbz;

  Vect coord, ncoord, ua, ue;
  double sangle;
  for(int i = 0; i < natom; i ++){
    coord = atom_init[i].coord;
    //apply pbc
    Vect dr = VSub(coord, cdisl[loopid]);
    PBCShift(coord, dr);
    //plastic displacement
    sangle = AtomSolidAngle(loopid, coord);
    ua     = VMul(b, -sangle/4./PI);
		//elastic displacement
		if(nu != 0){
			ue = ElasticDispLoop(coord, loopid);
		  ua = VAdd(ua, ue);
		}
    atom[i].coord = VAdd(atom[i].coord, ua);
    ApplyPBC(atom[i]);
    //to save the total displacement 
    u[i] = VAdd(u[i], ua); 
    //u[i] = VAdd(u[i], ue); 
  }
}
/*------------------------------------------------------*/
//to calculate the elastic displacement of a loop
/*------------------------------------------------------*/
Vect CreateDisl::ElasticDispLoop(Vect vP, int loopid)
{
  Vect vdisl01 = VSub(vdisl[loopid][1], vdisl[loopid][0]);
  Vect vdisl02 = VSub(vdisl[loopid][2], vdisl[loopid][0]);
  Vect cross   = VCross(vdisl01, vdisl02);
  Vect vdisl0a = VSub(vP, vdisl[loopid][0]);
  //if the vector 'cross' is not normalized, it may produce a vector with 
  //a large magnitude and affect the value of variable 'val' 
  cross   = VNorm(cross);
  vdisl0a = VNorm(vdisl0a);
  double val = VDot(cross, vdisl0a);

  Vect u, useg;
	InitV(u);
	Vect vA, vB, b;
  b = VMul(bv[loopid].v, bs[loopid]);

	for(int i = 0; i < ndisl[loopid]; i ++){
		vA = vdisl[loopid][i];
		if(i == ndisl[loopid] - 1) vB = vdisl[loopid][0];
		else vB = vdisl[loopid][i+1]; 
		useg = DispSegment(vP, vA, vB, b);
		u = VAdd(u, useg);
	}
	if(val > 0) return u;
	else return VMul(u, -1.);
}
/*------------------------------------------------------*/
Vect CreateDisl::ElasticDispTriangle(Vect vP, Vect vA, Vect vB, Vect vC, Vect b)
{
	Vect vAB = VSub(vB, vA);
	Vect vAC = VSub(vC, vA);
	Vect vAP = VSub(vP, vA);
	Vect cross = VCross(vAB, vAC);
	cross = VNorm(cross);
	double val = VDot(cross, VNorm(vAP));

	Vect u, uAB, uBC, uCA;
	InitV(u);
	uAB = DispSegment(vP, vA, vB, b);
	uBC = DispSegment(vP, vB, vC, b);
	uCA = DispSegment(vP, vC, vA, b);
	u = VAdd(u, uAB);
	u = VAdd(u, uBC);
	u = VAdd(u, uCA);
	if(val > 0) return u;
	else return VMul(u, -1.);
}
/*------------------------------------------------------*/
Vect CreateDisl::ElasticDispHeli(int helid, Vect coord, Vect b)
{
	Vect ue, ut;
	InitV(ue);
	for(int j = 0; j < ngridh[helid]; j ++){
		Vect vP = coord;
		Vect vA = gridh[helid][j].v0;
		Vect vB = gridh[helid][j].v1;
		Vect vC = gridh[helid][j].v2;
		ut = ElasticDispTriangle(vP, vA, vB, vC, b);
		ue = VAdd(ue, ut);
	}
	return ue;
}
/*------------------------------------------------------*/
void CreateDisl::ShiftAtoms()
{
	//to avoid atoms perfectly align with dislocation core
	Vect dr = {0.1341321, 0.19123481, 0.102348671};
	for(int i = 0; i < natom; i ++){
		atom[i].coord = VAdd(atom[i].coord, dr);
		ApplyPBC(atom[i]);
	}
}
/*------------------------------------------------------*/
Vect CreateDisl::DispSegment(Vect vP, Vect vA, Vect vB, Vect b)
{
	//vp: field point; vA, vB: end points of segment
	Vect vPA, vPB, tAB, nAB, u;
	vPA = VSub(vA, vP);
	vPB = VSub(vB, vP);

	tAB = VSub(vB, vA);
	tAB = VNorm(tAB);

	nAB = VCross(vPA, vPB);
	nAB = VNorm(nAB);

	Vect logterm;
	Vect angterm;
	logterm = VMul( VCross(b,tAB), log( (VLen(vPB)+VDot(vPB, tAB)) / (VLen(vPA)+VDot(vPA,tAB)) ) );
	logterm = VMul( logterm, -(1.-2.*nu) );
	Vect dv;
	dv      = VSub( VNorm(vPB), VNorm(vPA) );
	angterm = VMul( VCross(dv, nAB), VDot(b, nAB) );
	u = VAdd(logterm, angterm);
	u = VMul(u, 1./(8.*PI*(1.-nu)) );
	return u;
}
/*------------------------------------------------------*/
void CreateDisl::DiscretizeLoopTe(int tesfid)
{
  int nv = 3; //vertex number of each face     
  double dtheta = 2*PI / nv;
  //three vertices in local coordinate system
  for(int i = 0; i < nv; i ++){
    vtesf[tesfid][i].x = rtesf[tesfid] * cos(dtheta*i - PI/6.);
    vtesf[tesfid][i].y = rtesf[tesfid] * sin(dtheta*i - PI/6.);
    vtesf[tesfid][i].z = 0;
  }   

  Vect xn, yn, zn; //new axises 
  zn = normtesf[tesfid].v;
  yn =      lyt[tesfid].v;
  xn =      lxt[tesfid].v;
  //tranform into global coordinate system 
  double m[3][3], np[3], op[3];
  m[0][0] = VDot(xn, x.v)/VLen(xn)/VLen(x.v);
  m[1][0] = VDot(xn, y.v)/VLen(xn)/VLen(y.v);
  m[2][0] = VDot(xn, z.v)/VLen(xn)/VLen(z.v);

  m[0][1] = VDot(yn, x.v)/VLen(yn)/VLen(x.v);
  m[1][1] = VDot(yn, y.v)/VLen(yn)/VLen(y.v);
  m[2][1] = VDot(yn, z.v)/VLen(yn)/VLen(z.v);

  m[0][2] = VDot(zn, x.v)/VLen(zn)/VLen(x.v);
  m[1][2] = VDot(zn, y.v)/VLen(zn)/VLen(y.v);
  m[2][2] = VDot(zn, z.v)/VLen(zn)/VLen(z.v);

  for(int n = 0; n < nv; n ++){
    for(int i = 0; i < 3; i ++) np[i] = 0; 
    op[0] = vtesf[tesfid][n].x; 
    op[1] = vtesf[tesfid][n].y; 
    op[2] = vtesf[tesfid][n].z;
    for(int i = 0; i < 3; i ++){
      for(int j = 0; j < 3; j ++){
        np[i] += m[i][j] * op[j];
      }
    }
    vtesf[tesfid][n].x = np[0]; 
    vtesf[tesfid][n].y = np[1]; 
    vtesf[tesfid][n].z = np[2];
    vtesf[tesfid][n] = VAdd(vtesf[tesfid][n], ctesf[tesfid]);
  }

  //coordination of the forth vertex of the tetrahedron
  Vect center;
  center.x = (vtesf[tesfid][0].x + vtesf[tesfid][1].x + vtesf[tesfid][2].x)/3.;
  center.y = (vtesf[tesfid][0].y + vtesf[tesfid][1].y + vtesf[tesfid][2].y)/3.;
  center.z = (vtesf[tesfid][0].z + vtesf[tesfid][1].z + vtesf[tesfid][2].z)/3.;

  double a = sqrt(3.)*rtesf[tesfid];
  double h = sqrt(2.) / sqrt(3.) * a;
  vtesf[tesfid][3].x  = center.x - h * m[0][2];
  vtesf[tesfid][3].y  = center.y - h * m[1][2];
  vtesf[tesfid][3].z  = center.z - h * m[2][2];
}
/*------------------------------------------------------*/
double CreateDisl::AtomSolidAngleTe(int tesfid, int fid, Vect coord)
{
  int vface[4][3] = {{0, 1, 2},{0, 2, 3},{0, 1, 3},{1, 2, 3}}; //vertex index of each face
  int vertex[4] = {3, 1, 2, 0};                                //the forth vertex of tesf
  //check if the atom is on the dislocaiton plane  
  Vect vtesf01, vtesf02, vtesf03, vtesf0a, cross;
  int v0, v1, v2, v3;
  v0 = vface[fid][0];
  v1 = vface[fid][1];
  v2 = vface[fid][2];
  v3 = vertex[fid]; 

  vtesf01 = VSub(vtesf[tesfid][v1], vtesf[tesfid][v0]);
  vtesf02 = VSub(vtesf[tesfid][v2], vtesf[tesfid][v0]);
  vtesf03 = VSub(vtesf[tesfid][v3], vtesf[tesfid][v0]);
  vtesf0a = VSub(coord, vtesf[tesfid][v0]);
  cross   = VCross(vtesf01, vtesf02);
  cross   = VNorm(cross);
  vtesf0a = VNorm(vtesf0a);
  vtesf03 = VNorm(vtesf03);
  double val_a =  VDot(cross, vtesf0a);
  double val_v3 = VDot(cross, vtesf03);

  //case 1: atom on the dislocation plane
  Tetrahedron te;
  Triangle t;
  double sAngle = 0;
  if(fabs(val_a) < 1e-6){ 
    t.v0 = vtesf[tesfid][v0];
    t.v1 = vtesf[tesfid][v1];
    t.v2 = vtesf[tesfid][v2];
    sAngle = SolidAngleTriangle(t, coord); 
    return sAngle;
  }
  //case 2: atom lays out of the dislocation plane
  else{
    te.v0 = coord; 
    te.v1 = vtesf[tesfid][v0];
    te.v2 = vtesf[tesfid][v1];
    te.v3 = vtesf[tesfid][v2];
    sAngle = SolidAngleTetrahedron(te);
    if(val_a * val_v3 < 0) return sAngle; //atom and the 4th vertex are in the same side of the face
    else return -sAngle;
  }
}
/*------------------------------------------------------*/
void CreateDisl::ApplyDisplacementTe(int tesfid, int fid)
{
  //sequences of the vertex and the outward normls of each face follow the right hand rule
  //int vface[4][3] = {{0, 1, 2},{0, 2, 3},{0, 3, 1},{1, 3, 2}};
  int vface[4][3] = {{0, 1, 2},{0, 2, 3},{0, 1, 3},{1, 2, 3}}; //vertex index of each face
  int vertex[4] = {3, 1, 2, 0};                                // the left vertex

  int v0, v1, v2, v3;
  v0 = vface[fid][0];
  v1 = vface[fid][1];
  v2 = vface[fid][2];
  v3 = vertex[fid];

  Vect norm, vtesf01, vtesf02, vtesf03; 
  vtesf01 = VSub(vtesf[tesfid][v1], vtesf[tesfid][v0]);
  vtesf02 = VSub(vtesf[tesfid][v2], vtesf[tesfid][v0]);
  vtesf03 = VSub(vtesf[tesfid][v3], vtesf[tesfid][v0]);
  norm    = VCross(vtesf01, vtesf02);
  double vdot = VDot(norm, vtesf03);
  // normal of the face is NOT outward.
  if(vdot > 0) norm = VMul(norm, -1.); 
  NormalizeV(norm);

  Vect b;
  if(fid == 0){      // slip vector for Frank loop: <111>/3
    b = VMul(norm, sqrt(3.)/3.*lat_a); 
  }
  else if(fid > 0){  // slip vector for Shockely partial loop: <112>/6
    Vect mid = VMid(vtesf[tesfid][v0], vtesf[tesfid][v1]);
    Vect dir = VSub(vtesf[tesfid][v2], mid);
    NormalizeV(dir);
    b = VMul(dir, -sqrt(6.)/6.*lat_a);
  }

  Vect coord, ua;
  double sangle;
  for(int i = 0; i < natom; i ++){
    coord = atom_init[i].coord;
    //apply pbc 
    Vect dr = VSub(coord, ctesf[tesfid]);
    PBCShift(coord, dr);
		//plastic displacement
    sangle = AtomSolidAngleTe(tesfid, fid, coord);
    ua     = VMul(b, -sangle/4./PI);
		//elastic displacement
		if(nu != 0){
			Vect vP = coord;
			Vect vA = vtesf[tesfid][v0];
			Vect vB = vtesf[tesfid][v1];
			Vect vC = vtesf[tesfid][v2];
			Vect ue = ElasticDispTriangle(vP,vA,vB,vC,b);
			ua = VAdd(ua, ue);
		}
    atom[i].coord = VAdd(atom[i].coord, ua);
    ApplyPBC(atom[i]);
		//save displacement
    u[i]   = VAdd(u[i], ua);
  }
}
/*------------------------------------------------------*/
void CreateDisl::ReadGrid(int helid)
{
  char *line;
  line = new char[1024];
  FILE *fp; 
  char *first, *next, *ptr;
  int nitem, narg; 
  int rowNum = 0;
  char **items = new char* [32*sizeof(char*)]; // assume that 20 items at most each line

  fp = fopen(fheli[helid], "r");
  if(fp == NULL){
    printf("CanNOT find input file: %s\n", fheli[helid]);
    exit(1);
  }
  // default: there is no dislocation or stacking fault tetrahedron
  while(1) {  
    nitem = 0;
    if(fgets(line, 1024, fp) == NULL) { break; }
    first = nextword(line, &next);
    if(first == NULL) continue;           // blank line  
    rowNum ++;                            // add row number by one 
    ptr = next;
    items[0] = first;
    while(ptr) {
      nitem ++;
      items[nitem] = nextword(ptr, &next);
      if(!items[nitem]) break;
      ptr = next;
    }    
    narg = nitem;

    if(rowNum == 1) {
      ngridh[helid] = atoi(items[0]); //number of triangulation mesh 
      gridh[helid]  = new Triangle [ngridh[helid]];
    }
    else if(rowNum > 1){
      int id = (rowNum - 2) / 3; //triangle id in mesh
      if((rowNum - 2) % 3 == 0) {
        gridh[helid][id].v0.x = atof(items[0]); 
        gridh[helid][id].v0.y = atof(items[1]); 
        gridh[helid][id].v0.z = atof(items[2]); 
      }
      else if((rowNum - 2)%3 == 1){
        gridh[helid][id].v1.x = atof(items[0]); 
        gridh[helid][id].v1.y = atof(items[1]); 
        gridh[helid][id].v1.z = atof(items[2]); 
      }
      else if((rowNum - 2)%3 == 2){
        gridh[helid][id].v2.x = atof(items[0]); 
        gridh[helid][id].v2.y = atof(items[1]); 
        gridh[helid][id].v2.z = atof(items[2]); 
      }
    }
  }
  fclose(fp);
}
/*------------------------------------------------------*/
double CreateDisl::AtomSolidAngleHeli(int helid, Vect coord)
{
  double sAngle = 0;
  Vect norm;
  Vect v0a;
  Vect v01, v02;
  double val;
  Tetrahedron te; 

  for(int j = 0; j < ngridh[helid]; j ++){
    v01  = VSub(gridh[helid][j].v1, gridh[helid][j].v0);
    v02  = VSub(gridh[helid][j].v2, gridh[helid][j].v0);
    v0a  = VSub(coord,              gridh[helid][j].v0);
    norm = VCross(v01, v02);
    norm = VNorm(norm);
    v0a  = VNorm(v0a);
    double val  = VDot(v0a, norm);
    //atom on the plane
    if(fabs(val) < 1e-6){
      sAngle += SolidAngleTriangle(gridh[helid][j], coord);
    }
    else{
      te.v0 = coord;
      te.v1 = gridh[helid][j].v0;
      te.v2 = gridh[helid][j].v1;
      te.v3 = gridh[helid][j].v2;
      if(val > 0)
        sAngle +=  SolidAngleTetrahedron(te);
      else 
        sAngle += -SolidAngleTetrahedron(te);
    }
  }
  return sAngle;
}
/*------------------------------------------------------*/
void CreateDisl::ApplyDisplacementHeli(int helid)
{
  Vect b;
  b.x   = bvh[helid].v.x * bsh[helid];
  b.y   = bvh[helid].v.y * bsh[helid];
  b.z   = bvh[helid].v.z * bsh[helid];
  //transform Burgers vector to Cartesian coordinate system
  double vlen = VLen(b);
  double cosbx = VDot(b, x.v)/vlen/VLen(x.v);
  double cosby = VDot(b, y.v)/vlen/VLen(y.v);
  double cosbz = VDot(b, z.v)/vlen/VLen(z.v);
  b.x = vlen * cosbx;
  b.y = vlen * cosby;
  b.z = vlen * cosbz;

  Vect coord, ua, ue;
  double sangle;
  for(int i = 0; i < natom; i ++){
    coord = atom_init[i].coord;
		//plastic displacement
    sangle = AtomSolidAngleHeli(helid, coord);
    ua = VMul(b, -sangle/4./PI);
		//elastic displacement
		if(nu != 0){
			ue = ElasticDispHeli(helid, coord, b);
			ua = VAdd(ua, ue);
		}
    atom[i].coord = VAdd(atom[i].coord, ua);
    ApplyPBC(atom[i]);
		//save total displacement
		u[i] = VAdd(u[i], ua);
  }
}
/*------------------------------------------------------*/
/*------------------------------------------------------*/
void CreateDisl::Output(char *file)
{
  //reset box size slightly when non-pbc is applied 
  if(pbcx != 1 && resetpbcx != 1) {
    regInfo.xlo += -0.5; 
    regInfo.xhi +=  0.5;
    regInfo.xprd = regInfo.xhi - regInfo.xlo;
  }
  if(pbcy != 1 && resetpbcy != 1) {
    regInfo.ylo += -0.5; 
    regInfo.yhi +=  0.5;
    regInfo.yprd = regInfo.yhi - regInfo.ylo;
  }
  if(pbcz != 1 && resetpbcz != 1) {
    regInfo.zlo += -0.5; 
    regInfo.zhi +=  0.5;
    regInfo.zprd = regInfo.zhi - regInfo.zlo;
  }

  // output in lammps data-file format
  FILE *fp;
  char *fname = new char [64];
  strcpy(fname, file);
  strcat(fname, ".lmp");
  fp = fopen(fname, "w");
  fprintf(fp, "LAMMPS description\n");
  fprintf(fp, "\n");
  fprintf(fp, "%d atoms\n", natom - natom_del);
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
  for(int i = 0; i < natom; i ++){
    if(atom[i].stat == 0) continue;
    fprintf(fp, "%d %d %f %f %f\n", id, atom[i].type, 
        atom[i].coord.x, atom[i].coord.y, atom[i].coord.z);
    id ++;
  }
  fclose(fp);

  printf("---------------Output file------------------\n");
  printf("  %s: in lammps data-file format\n", fname);
  printf("--------------------------------------------\n");

  //dump_file
  strcpy(fname, file);
  strcat(fname,".dump");
  fp = fopen(fname, "w");
  fprintf(fp, "ITEM: TIMESTEP\n");
  fprintf(fp, "0\n");
  fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
  fprintf(fp, "%d\n", natom - natom_del);
  //fprintf(fp, "%d\n", natom - natom_del + ngridh[0]);
  fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n");
  fprintf(fp, "%f %f\n", regInfo.xlo, regInfo.xhi);
  fprintf(fp, "%f %f\n", regInfo.ylo, regInfo.yhi);
  fprintf(fp, "%f %f\n", regInfo.zlo, regInfo.zhi);
  fprintf(fp, "ITEM: ATOMS id type x y z ux uy uz\n");
  id = 1;
  for(int i = 0; i < natom; i ++){
    if(atom[i].stat == 0) continue;
    fprintf(fp, "%d %d %f %f %f %f %f %f\n", id, atom[i].type,
        atom[i].coord.x, atom[i].coord.y, atom[i].coord.z, u[i].x, u[i].y, u[i].z);
    id ++;
  }
  fclose(fp);

  printf("---------------Output file------------------\n");
  printf("  %s: in lammps dump-file format\n", fname);
  printf("--------------------------------------------\n");

	/*
	//===================================================================
	//inital atom pos
  strcpy(fname, file);
  strcat(fname, ".lmp.init");
  fp = fopen(fname, "w");
  fprintf(fp, "LAMMPS description\n");
  fprintf(fp, "\n");
  fprintf(fp, "%d atoms\n", natom - natom_del);
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
  id = 1;
  for(int i = 0; i < natom; i ++){
    if(atom[i].stat == 0) continue;
    fprintf(fp, "%d %d %f %f %f\n", id, atom_init[i].type, 
        atom_init[i].coord.x, atom_init[i].coord.y, atom_init[i].coord.z);
    id ++;
  }
  fclose(fp);

  printf("---------------Output file------------------\n");
  printf("  %s: in lammps data-file format (initial position)\n", fname);
  printf("--------------------------------------------\n");
	//===================================================================
	*/
}
/*------------------------------------------------------*/
void CreateDisl::ReadConfig(char *file)
{
  char *line, *copy;
  line = new char[1024];
  FILE *fp; 
  char *first, *next, *ptr;
  int nitem, narg, rowNum = 0;
  char **items = new char* [32*sizeof(char*)]; // assume that 32 items at most each line

  ntype = 0;    //used to store the number of atom types
  fp = fopen(file, "r");
  if(fp == NULL){
    printf("CanNOT find input file: %s\n", file);
    exit(1);
  }
  while(1) {  
    nitem = 0;
    if(fgets(line, 1024, fp) == NULL) break;
    first = nextword(line, &next);
    if(first == NULL) continue; // blank line  
    rowNum ++;                  // add row number by one 
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
      u = new Vect[natom];
      for(int i = 0; i < natom; i ++){
        u[i].x = 0;
        u[i].y = 0;
        u[i].z = 0;
      }
    }
    if(rowNum == 3){
      ntype = atoi(items[0]);
      // reset type of atoms around dislocation
      // ntype = ntype;
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

	//to avoid atoms perfectly align with dislocaiton core 
	//when calculat elastic displacment
	if(atomShiftFlag == 1) ShiftAtoms();

  //store the initial position of atoms
  atom_init = new Atom [natom];
  for(int i = 0; i < natom; i ++){
    atom_init[i] = atom[i];
  }
}
/*------------------------------------------------------*/
char* CreateDisl::nextword(char *str, char **next)
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
/*-----------------------------------------------------------------------
  Delete the atoms that lay close to the neighbors after applying displacment 
-------------------------------------------------------------------------*/
void CreateDisl::DeleteAtom(int loopid)
{
  // set flags for each atom to indicate whether it is deleted or NOT
  rCut = 2.0; // if the distance between two atoms is less than rCut, randomly delete one of them
  nebr = NULL;
  cellList = NULL;
  // subdivide the region into cells
  VectI cc, m1v, m2v, cells;
  Vect shift;
  Vect dr, rs;
  int N_OFFSET = 14;
  VectI vOff[] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{-1,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},{-1,1,1},
        {-1,0,1},{-1,-1,1},{0,-1,1},{1,-1,1}};
  double cellLength = rCut * 1.5;
  cells.x = int(regInfo.xprd / cellLength);
  cells.y = int(regInfo.yprd / cellLength);
  cells.z = int(regInfo.zprd / cellLength);

  cellListLen = natom + cells.x * cells.y * cells.z;
  cellList = new int[cellListLen];
  for(int i = natom; i < cellListLen; i ++) cellList[i] = -1;
  for(int i = 0; i < natom; i ++) {
    rs.x = atom[i].coord.x + regInfo.xprd * 0.5;
    rs.y = atom[i].coord.y + regInfo.yprd * 0.5;
    rs.z = atom[i].coord.z + regInfo.zprd * 0.5;
    cc.x = int (rs.x / cellLength);
    cc.y = int (rs.y / cellLength);
    cc.z = int (rs.z / cellLength);
    if(cc.x >= cells.x) cc.x = cells.x - 1;
    if(cc.x <= 0) cc.x = 0;
    if(cc.y >= cells.y) cc.y = cells.y - 1;
    if(cc.y <= 0) cc.y = 0;
    if(cc.z >= cells.z) cc.z = cells.z - 1;
    if(cc.z <= 0) cc.z = 0;
    int c = cc.z * cells.x * cells.y + cc.y * cells.x + cc.x + natom;
    if(c > cellListLen) { printf("BUILDING CELLIST ERROR!\n"); return;}
    cellList[i] = cellList[c];
    cellList[c] = i;
  }

  // build neighbor 
  nebrListLen = 0;
  Nebr *tmp = NULL;
  //================================
  for(int m1x = 0; m1x < cells.x; m1x ++) {
    for(int m1y = 0; m1y < cells.y; m1y ++) {
      for(int m1z = 0; m1z < cells.z; m1z ++) {
        m1v.x = m1x;
        m1v.y = m1y;
        m1v.z = m1z;
        int m1 = m1v.z * cells.x * cells.y + m1v.y * cells.x + m1v.x+ natom;
        for(int offset = 0; offset < N_OFFSET; offset ++) {
          m2v.x = m1v.x + vOff[offset].x;
          m2v.y = m1v.y + vOff[offset].y;
          m2v.z = m1v.z + vOff[offset].z;
          shift.x = 0; shift.y = 0; shift.z = 0;
          //-------apply periodic boundary condition-----        
          if(m2v.x >= cells.x){ m2v.x = 0; shift.x = regInfo.xprd; } 
          else if(m2v.x < 0){ m2v.x = cells.x - 1; shift.x = -regInfo.xprd; }
          if(m2v.y >= cells.y){ m2v.y = 0; shift.y = regInfo.yprd; }
          else if(m2v.y < 0){ m2v.y = cells.y - 1; shift.y = -regInfo.yprd; }
          if(m2v.z >= cells.z){ m2v.z = 0; shift.z = regInfo.zprd; }
          else if(m2v.z < 0){ m2v.z = cells.z - 1; shift.z = -regInfo.zprd;
          }
          //--------end of periodic boundary condition-------
          int m2 = m2v.z * cells.x * cells.y + m2v.y * cells.x + m2v.x + natom;
          for(int j1 = cellList[m1]; j1 >= 0; j1 = cellList[j1]) {
            for(int j2 = cellList[m2]; j2 >= 0; j2 = cellList[j2]) {
              if(m1 != m2 || j2 < j1) { dr.x = atom[j1].coord.x - atom[j2].coord.x;
                dr.y = atom[j1].coord.y - atom[j2].coord.y;
                dr.z = atom[j1].coord.z - atom[j2].coord.z;
              
                dr.x -= shift.x;
                dr.y -= shift.y;
                dr.z -= shift.z;              

                double rr = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
                if(rr < rCut * rCut) {
                  nebrListLen ++;
                  if(nebrListLen >= natom) { 
                    printf("Array for nebrlist NOT enough\n");
                    exit(1);
                  } 
                  tmp = (Nebr*)realloc(nebr, nebrListLen*sizeof(Nebr));  
                  if(tmp == NULL){
                    printf("Memory allocaton for nebr-list error!\n"); exit(1);
                  }
                  nebr = tmp;
                  nebr[nebrListLen - 1].i = j1;
                  nebr[nebrListLen - 1].j = j2;
                }
              }// end of if          
            } // end of for
          } // end of for
        }// end of for
      }// end of for
    }
  }// end of for
  //================================
  // deleting operation begions
  //================================
  int ni, nj;  
  double proi, proj; 
  int ndel = 0;
  //natom_del = 0;
  Vect vi, vj, norm, ndisl;
  for(int i = 0; i < nebrListLen; i ++) {
    ni = nebr[i].i; 
    nj = nebr[i].j;
    if(atom[ni].stat == 1 && atom[nj].stat == 1) {
      //to distinguish atom planes on the reference configuraiton
      vi = atom_init[ni].coord;
      vj = atom_init[nj].coord;
      //apply pbc
      Vect dri = VSub(vi, cdisl[loopid]);
      Vect drj = VSub(vj, cdisl[loopid]);
      PBCShift(vi, dri);
      PBCShift(vj, drj);
      
      ndisl = normdisl[loopid].v;
      norm.x = VDot(ndisl, x.v)/VLen(ndisl)/VLen(x.v);
      norm.y = VDot(ndisl, y.v)/VLen(ndisl)/VLen(y.v);
      norm.z = VDot(ndisl, z.v)/VLen(ndisl)/VLen(z.v);

      proi = VDot(vi, norm);
      proj = VDot(vj, norm);
      if(proi >= proj)  atom[nj].stat = 0; 
      else if(proi < proj)  atom[ni].stat = 0; 
      ndel ++;
    }
  }
  natom_del += ndel; // add to the total number
  printf("Number of deleted atoms: %d\n", ndel);
  //release memory
  if(nebr != NULL) free (nebr);
  if(cellList != NULL) free(cellList);
}
/*-----------------------------------------------------------------------
  Delete atoms that lay close to the neighbors after applying displacment 
-------------------------------------------------------------------------*/
void CreateDisl::DeleteAtomTe(int tesfid, int fid)
{
  // set flags for each atom to indicate whether it is deleted or NOT
  rCut = 1.7; // if distance between two atoms less than rCut, randomly delete one of them
  nebr = NULL;
  cellList = NULL;
  // cell subdivide the region
  VectI cc, m1v, m2v, cells;
  Vect shift;
  Vect dr, rs;
  int N_OFFSET = 14;
  VectI vOff[] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{-1,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},{-1,1,1},
        {-1,0,1},{-1,-1,1},{0,-1,1},{1,-1,1}};
  double cellLength = rCut * 1.5;
  cells.x = int(regInfo.xprd / cellLength);
  cells.y = int(regInfo.yprd / cellLength);
  cells.z = int(regInfo.zprd / cellLength);

  cellListLen = natom + cells.x * cells.y * cells.z;
  cellList = new int[cellListLen];

  for(int i = natom; i < cellListLen; i ++) cellList[i] = -1;
  for(int i = 0; i < natom; i ++) {
    rs.x = atom[i].coord.x + regInfo.xprd * 0.5;
    rs.y = atom[i].coord.y + regInfo.yprd * 0.5;
    rs.z = atom[i].coord.z + regInfo.zprd * 0.5;
    cc.x = int (rs.x / cellLength);
    cc.y = int (rs.y / cellLength);
    cc.z = int (rs.z / cellLength);
    if(cc.x >= cells.x) cc.x = cells.x - 1;
    if(cc.x <= 0) cc.x = 0;
    if(cc.y >= cells.y) cc.y = cells.y - 1;
    if(cc.y <= 0) cc.y = 0;
    if(cc.z >= cells.z) cc.z = cells.z - 1;
    if(cc.z <= 0) cc.z = 0;
    int c = cc.z * cells.x * cells.y + cc.y * cells.x + cc.x + natom;
    if(c > cellListLen) { printf("BUILDING CELLIST ERROR!\n"); return;}
    cellList[i] = cellList[c];
    cellList[c] = i;
  }

  // build neighbor 
  nebrListLen = 0;
  Nebr *tmp = NULL;
  //================================
  for(int m1x = 0; m1x < cells.x; m1x ++) {
    for(int m1y = 0; m1y < cells.y; m1y ++) {
      for(int m1z = 0; m1z < cells.z; m1z ++) {
        m1v.x = m1x;
        m1v.y = m1y;
        m1v.z = m1z;
        int m1 = m1v.z * cells.x * cells.y + m1v.y * cells.x + m1v.x+ natom;
        for(int offset = 0; offset < N_OFFSET; offset ++) {
          m2v.x = m1v.x + vOff[offset].x;
          m2v.y = m1v.y + vOff[offset].y;
          m2v.z = m1v.z + vOff[offset].z;
          shift.x = 0; shift.y = 0; shift.z = 0;
          //-------apply periodic boundary condition-----        
          if(m2v.x >= cells.x){ m2v.x = 0; shift.x = regInfo.xprd; } 
          else if(m2v.x < 0){ m2v.x = cells.x - 1; shift.x = -regInfo.xprd; }
          if(m2v.y >= cells.y){ m2v.y = 0; shift.y = regInfo.yprd; }
          else if(m2v.y < 0){ m2v.y = cells.y - 1; shift.y = -regInfo.yprd; }
          if(m2v.z >= cells.z){ m2v.z = 0; shift.z = regInfo.zprd; }
          else if(m2v.z < 0){ m2v.z = cells.z - 1; shift.z = -regInfo.zprd;
          }
          //--------end of periodic boundary condition-------
          int m2 = m2v.z * cells.x * cells.y + m2v.y * cells.x + m2v.x + natom;
          for(int j1 = cellList[m1]; j1 >= 0; j1 = cellList[j1]) {
            for(int j2 = cellList[m2]; j2 >= 0; j2 = cellList[j2]) {
              if(m1 != m2 || j2 < j1) { 
                dr.x = atom[j1].coord.x - atom[j2].coord.x;
                dr.y = atom[j1].coord.y - atom[j2].coord.y;
                dr.z = atom[j1].coord.z - atom[j2].coord.z;
                dr.x -= shift.x;
                dr.y -= shift.y;
                dr.z -= shift.z;              

                double rr = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
                if(rr < rCut * rCut) {
                  nebrListLen ++;
                  if(nebrListLen >= natom) { 
                    printf("Array for nebrlist NOT enough\n");
                    exit(1);
                  } 
                  tmp = (Nebr*)realloc(nebr, nebrListLen*sizeof(Nebr));  
                  if(tmp == NULL){
                    printf("Memory allocaton for nebr-list error!\n"); exit(1);
                  }
                  nebr = tmp;
                  nebr[nebrListLen - 1].i = j1;
                  nebr[nebrListLen - 1].j = j2;
                }
              }// end of if          
            } // end of for
          } // end of for
        }// end of for
      }// end of for
    }
  }// end of for
  //================================
  // deleting operation begions
  //================================
  //normal
  int vface[4][3] = {{0, 1, 2},{0, 2, 3},{0, 3, 1},{1, 3, 2}};
  int v0, v1, v2;

  Vect norm, vtesf01, vtesf02; 
  v0 = vface[fid][0];
  v1 = vface[fid][1];
  v2 = vface[fid][2];
  vtesf01 = VSub(vtesf[tesfid][v1], vtesf[tesfid][v0]);
  vtesf02 = VSub(vtesf[tesfid][v2], vtesf[tesfid][v0]);
  norm    = VCross(vtesf01, vtesf02);
  NormalizeV(norm);

  int ni, nj;  
  double proi, proj; 
  int ndel = 0;
  //natom_del = 0;
  Vect vi, vj;
  for(int i = 0; i < nebrListLen; i ++) {
    ni = nebr[i].i; 
    nj = nebr[i].j;
    if(atom[ni].stat == 1 && atom[nj].stat == 1) {
      //to distinguish atom planes on the reference configuraiton
      vi = atom_init[ni].coord;
      vj = atom_init[nj].coord;
      //apply pbc
      Vect dri = VSub(vi, ctesf[tesfid]);
      Vect drj = VSub(vj, ctesf[tesfid]);
      PBCShift(vi, dri);
      PBCShift(vj, drj);

      proi = VDot(vi, norm);
      proj = VDot(vj, norm);
      if(proi >= proj) {
        atom[nj].stat = 0;
        ndel ++;
      }
      else if(proi < proj) { 
        atom[ni].stat = 0;
        ndel ++;
      }
    }
  }
  natom_del += ndel; // add to the total number
  printf("Number of deleted atoms: %d\n", ndel);
  //release memory
  if(nebr != NULL) free(nebr);
  if(cellList != NULL) free(cellList);
}
/*-----------------------------------------------------------------------
  Delete atoms that lay close to the neighbors after applying displacment 
-------------------------------------------------------------------------*/
void CreateDisl::DeleteAtomHeli(int helid)
{
  // set flags for each atom to indicate whether it is deleted or NOT
  rCut = 2.0; // if distance between two atoms less than rCut, randomly delete one of them
  nebr = NULL;
  cellList = NULL;
  // cell subdivide the region
  VectI cc, m1v, m2v, cells;
  Vect shift;
  Vect dr, rs;
  int N_OFFSET = 14;
  VectI vOff[] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{-1,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},{-1,1,1},
        {-1,0,1},{-1,-1,1},{0,-1,1},{1,-1,1}};
  double cellLength = rCut * 1.5;
  cells.x = int(regInfo.xprd / cellLength);
  cells.y = int(regInfo.yprd / cellLength);
  cells.z = int(regInfo.zprd / cellLength);

  cellListLen = natom + cells.x * cells.y * cells.z;
  cellList = new int[cellListLen];
  for(int i = natom; i < cellListLen; i ++) cellList[i] = -1;
  for(int i = 0; i < natom; i ++) {
    rs.x = atom[i].coord.x + regInfo.xprd * 0.5;
    rs.y = atom[i].coord.y + regInfo.yprd * 0.5;
    rs.z = atom[i].coord.z + regInfo.zprd * 0.5;
    cc.x = int (rs.x / cellLength);
    cc.y = int (rs.y / cellLength);
    cc.z = int (rs.z / cellLength);
    if(cc.x >= cells.x) cc.x = cells.x - 1;
    if(cc.x <= 0) cc.x = 0;
    if(cc.y >= cells.y) cc.y = cells.y - 1;
    if(cc.y <= 0) cc.y = 0;
    if(cc.z >= cells.z) cc.z = cells.z - 1;
    if(cc.z <= 0) cc.z = 0;
    int c = cc.z * cells.x * cells.y + cc.y * cells.x + cc.x + natom;
    if(c > cellListLen) { printf("BUILDING CELLIST ERROR!\n"); return;}
    cellList[i] = cellList[c];
    cellList[c] = i;
  }

  // build neighbor 
  nebrListLen = 0;
  Nebr *tmp = NULL;
  //================================
  for(int m1x = 0; m1x < cells.x; m1x ++) {
    for(int m1y = 0; m1y < cells.y; m1y ++) {
      for(int m1z = 0; m1z < cells.z; m1z ++) {
        m1v.x = m1x;
        m1v.y = m1y;
        m1v.z = m1z;
        int m1 = m1v.z * cells.x * cells.y + m1v.y * cells.x + m1v.x+ natom;
        for(int offset = 0; offset < N_OFFSET; offset ++) {
          m2v.x = m1v.x + vOff[offset].x;
          m2v.y = m1v.y + vOff[offset].y;
          m2v.z = m1v.z + vOff[offset].z;
          shift.x = 0; shift.y = 0; shift.z = 0;
          //-------apply periodic boundary condition-----        
          if(m2v.x >= cells.x){ m2v.x = 0; shift.x = regInfo.xprd; } 
          else if(m2v.x < 0){ m2v.x = cells.x - 1; shift.x = -regInfo.xprd; }
          if(m2v.y >= cells.y){ m2v.y = 0; shift.y = regInfo.yprd; }
          else if(m2v.y < 0){ m2v.y = cells.y - 1; shift.y = -regInfo.yprd; }
          if(m2v.z >= cells.z){ m2v.z = 0; shift.z = regInfo.zprd; }
          else if(m2v.z < 0){ m2v.z = cells.z - 1; shift.z = -regInfo.zprd; }
          //--------end of periodic boundary condition-------
          int m2 = m2v.z * cells.x * cells.y + m2v.y * cells.x + m2v.x + natom;
          for(int j1 = cellList[m1]; j1 >= 0; j1 = cellList[j1]) {
            for(int j2 = cellList[m2]; j2 >= 0; j2 = cellList[j2]) {
              if(m1 != m2 || j2 < j1) { 
                dr.x = atom[j1].coord.x - atom[j2].coord.x;
                dr.y = atom[j1].coord.y - atom[j2].coord.y;
                dr.z = atom[j1].coord.z - atom[j2].coord.z;
              
                dr.x -= shift.x;
                dr.y -= shift.y;
                dr.z -= shift.z;              

                double rr = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
                if(rr < rCut * rCut) {
                  nebrListLen ++;
                  if(nebrListLen >= natom) { 
                    printf("Array for nebrlist NOT enough\n");
                    exit(1);
                  } 
                  tmp = (Nebr*)realloc(nebr, nebrListLen*sizeof(Nebr));  
                  if(tmp == NULL){
                    printf("Memory allocaton for nebr-list error!\n"); exit(1);
                  }
                  nebr = tmp;
                  nebr[nebrListLen - 1].i = j1;
                  nebr[nebrListLen - 1].j = j2;
                }
              }// end of if          
            } // end of for
          } // end of for
        }// end of for
      }// end of for
    }
  }// end of for
  //================================
  // deleting operation begions
  //================================
  int ni, nj;  
  double proi, proj; 
  int ndel = 0;
  //natom_del = 0;
  Vect vi, vj;
  Vect b; //direction of b in global coordinate system
  for(int i = 0; i < nebrListLen; i ++) {
    ni = nebr[i].i; 
    nj = nebr[i].j;
    if(atom[ni].stat == 1 && atom[nj].stat == 1) {
      //to distinguish atom planes on the reference configuraiton
      vi = atom_init[ni].coord;
      vj = atom_init[nj].coord;

      b.x = VDot(bvh[helid].v, x.v)/VLen(bvh[helid].v)/VLen(x.v);
      b.y = VDot(bvh[helid].v, y.v)/VLen(bvh[helid].v)/VLen(y.v);
      b.z = VDot(bvh[helid].v, z.v)/VLen(bvh[helid].v)/VLen(z.v);

      proi = VDot(vi, b);
      proj = VDot(vj, b);
      if(proi >= proj)  atom[nj].stat = 0; 
      else if(proi < proj)  atom[ni].stat = 0; 
      ndel ++;
    }
  }
  natom_del += ndel; // add to the total number
  printf("Number of deleted atoms: %d\n", ndel);
  //release memory
  if(nebr != NULL) free (nebr);
  if(cellList != NULL) free(cellList);
}
//-----------------------------------------------------------------------
// Add atoms for extrinsic Frank dislocaton
//-----------------------------------------------------------------------
void CreateDisl::AddAtom(int loopid)
{
  Vect norm; //normal of the slip plane 
  norm.x = VDot(normdisl[loopid].v, x.v)/VLen(normdisl[loopid].v)/VLen(x.v);
  norm.y = VDot(normdisl[loopid].v, y.v)/VLen(normdisl[loopid].v)/VLen(y.v);
  norm.z = VDot(normdisl[loopid].v, z.v)/VLen(normdisl[loopid].v)/VLen(z.v);

  Vect   center = cdisl[loopid];
  Vect   dv; // vector: center_of_plane -> atom[i]
  Atom   *atomAdded = NULL;
  int    nAtomAdded = 0;
  Atom   *tmp;
  double dot;
//double b = fabs(VLen(bv[loopid].v)*bs[loopid]*lat_a); //size of burger's vector
  double b = fabs(VLen(bv[loopid].v)*bs[loopid]);       //size of burger's vector
  double up1, lower1;      //nearest atomic planes above/below the slip plane
  double up2, lower2;
  double mid;              //middle position between up1 and lower1     
  for(int i = 0; i < natom; i ++){
    dv  = VSub(atom_init[i].coord, center);
    dot = VDot(dv, norm);  //projection to the normal direction
    if(fabs(dot) < 1e-6 ||
       (fabs(dot) > 1e-6 && fabs(dot - b) > 1e-6  && dot > 0 && dot < b)) {
      up1 = dot;
      up2 = dot + b;
      lower1 = up1 - b;
      mid = (up1 + lower1) / 2.;
      break;
    }
  }

  //fetch the first atom layer above up1
  double move = mid - up2;
  Vect um; 
  Vect ud;
  Vect coord;
  for(int i = 0; i < natom; i ++){
    coord = atom_init[i].coord;
    dv  = VSub(coord, center);
    dot = VDot(dv, norm); //projection to the normal direction

    if(dot < 0) continue; //only consider atoms in one side of the slip plane
    if(!IsInsideLoop(loopid, coord)) continue; //if atoms lay outside the region

    if(fabs(dot - b) < 1e-6 ||
       (fabs(dot - b) > 1e-6 && fabs(dot - 2.0*b) > 1e-6  && dot > b && dot < 2.0*b) ) {
      nAtomAdded ++;
      tmp = (Atom*)realloc(atomAdded, nAtomAdded*sizeof(Atom));
      if(tmp == NULL){
        printf("Memory allocaiton error n AddAtom(int)\n"); exit(1);
      }
      atomAdded = tmp;
      atomAdded[nAtomAdded - 1] = atom_init[i];
      um = VMul(norm, move); // up2 -> middle plane between up1 and lower1
      coord = VAdd(atom_init[i].coord, um);
      atomAdded[nAtomAdded - 1].coord = coord; 
    }
  }
  // add the new atoms to Array atom[]
  for(int i = 0; i < nAtomAdded; i ++){
    natom ++;
    tmp = (Atom*)realloc(atom, natom * sizeof(Atom));
    if(tmp == NULL){
      printf("Memory allocaiton error in AddAtom(int)\n"); exit(1);
    }
    atom = tmp;
    atom[natom - 1] = atomAdded[i];
    atom[natom - 1].type = 1;
  }
  if(atomAdded != NULL) free(atomAdded);

  // array recording atom displacement should also be adjusted 
  Vect *utmp = NULL;
  utmp = (Vect*)realloc(u, natom * sizeof(Vect));
  if(utmp == NULL){
    printf("Memory allocaiton error in AddAtom(int)\n"); exit(1);
  }
  u = utmp;
  for(int i = natom - nAtomAdded; i < natom; i ++){
    u[i].x = 0;  
    u[i].y = 0;  
    u[i].z = 0;  
  }
  //release memory
  printf("Added atoms: %d\n", nAtomAdded);
}
/*------------------------------------------------------*/
int CreateDisl::IsInsideLoop(int loopid, Vect coord)
{
  Vect ve, va, vc;
  Vect center, norm, cross;
  center = cdisl[loopid];
  norm.x = VDot(normdisl[loopid].v, x.v)/VLen(normdisl[loopid].v)/VLen(x.v);
  norm.y = VDot(normdisl[loopid].v, y.v)/VLen(normdisl[loopid].v)/VLen(y.v);
  norm.z = VDot(normdisl[loopid].v, z.v)/VLen(normdisl[loopid].v)/VLen(z.v);

  int next; 
  int v0, v1, vt;
  vt = ndisl[loopid];
  int flag = 1;
  for(int i = 0; i < vt; i++){
    v0 = i;
    v1 = v0 + 1;
    if(v1 == vt) v1 = 0;
    ve = VSub(vdisl[loopid][v1], vdisl[loopid][v0]);
    va = VSub(coord,  vdisl[loopid][v0]);
    vc = VSub(center, vdisl[loopid][v0]);
    cross = VCross(norm, ve);
    if(VDot(cross, va) * VDot(cross, vc) < 0) flag = 0; 
  }
  if(flag == 1) return 1;
  else return 0;
}
/*------------------------------------------------------*/
void CreateDisl::ApplyPBC(Atom &a)
{
  // '1': pbc, '0': non-pbc
  if(pbcx == 1){ 
    if(a.coord.x < regInfo.xlo) a.coord.x += regInfo.xprd;
    else if(a.coord.x >= regInfo.xhi) a.coord.x -= regInfo.xprd;
  }
  if(pbcy == 1){
    if(a.coord.y < regInfo.ylo) a.coord.y += regInfo.yprd;
    else if(a.coord.y >= regInfo.yhi) a.coord.y -= regInfo.yprd;
  }
  if(pbcz == 1){
    if(a.coord.z < regInfo.zlo) a.coord.z += regInfo.zprd;
    else if(a.coord.z >= regInfo.zhi) a.coord.z -= regInfo.zprd;
  }
}
/*------------------------------------------------------*/
void CreateDisl::PBCShift(Vect &v, Vect dr)
{
  if(pbcx == 1){
    if(fabs(dr.x) > regInfo.xprd * 0.5){
      if(dr.x > 0) v.x -= regInfo.xprd;
      else         v.x += regInfo.xprd;
    }
  }
  if(pbcy == 1){
    if(fabs(dr.y) > regInfo.yprd * 0.5){
      if(dr.y > 0) v.y -= regInfo.yprd;
      else         v.y += regInfo.yprd;
    }
  }
  if(pbcz == 1){
    if(fabs(dr.z) > regInfo.zprd * 0.5){
      if(dr.z > 0) v.z -= regInfo.zprd;
      else         v.z += regInfo.zprd;
    }
  }
}
/*------------------------------------------------------*/
double CreateDisl::Fraction2Decimal(char *para)
{
  int Length = 32;
  //numerator and denominator
  char *nmr = new char[Length]; 
  char *dnm = new char[Length];

  int n = strcspn(para, "/");
  if(n == strlen(para)){ 
    return atof(para); //when '/' is NOT found 
  } else {
    strncpy(nmr, para, n);
    strcpy(dnm, &para[n+1]);
    return atof(nmr)/atof(dnm); //when '/' is found
  }
}
/*------------------------------------------------------*/
int CreateDisl::Orthogonal(Vect l, Vect m, Vect n)
{
  if(fabs(VDot(l,m)) > 1e-10) return 0;
  if(fabs(VDot(l,n)) > 1e-10) return 0;
  if(fabs(VDot(m,n)) > 1e-10) return 0;
  /*
  if(l.x * m.x + l.y * m.y + l.z * m.z) return 0;
  if(l.x * n.x + l.y * n.y + l.z * n.z) return 0;
  if(m.x * n.x + m.y * n.y + m.z * n.z) return 0;
  */
  return 1;
}
/*------------------------------------------------------*/
int CreateDisl::Right_handed(Vect l, Vect m, Vect n)
{
  Vect c; 
  c.x = l.y * m.z - l.z * m.y;
  c.y = l.z * m.x - l.x * m.z;
  c.z = l.x * m.y - l.y * m.x;
  if(c.x*n.x + c.y*n.y + c.z*n.z < 0) return 0;
  return 1;
}
/*------------------------------------------------------*/
//set the center of bottom face of tetrahedron
//to the position of its nearest atom, and shift off position 
//slightly 
/*------------------------------------------------------*/
void CreateDisl::FindNearestAtom()
{
  rCut = 0.5*lat_a; 
  nebr = NULL;
  cellList = NULL;
  VectI cc, m1v, m2v, cells;
  Vect shift;
  Vect dr, rs;
  int N_OFFSET = 14;
  VectI vOff[] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{-1,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},{-1,1,1},
        {-1,0,1},{-1,-1,1},{0,-1,1},{1,-1,1}};
  double cellLength = rCut * 1.5;
  cells.x = int(regInfo.xprd / cellLength);
  cells.y = int(regInfo.yprd / cellLength);
  cells.z = int(regInfo.zprd / cellLength);

  cellListLen = natom + cells.x * cells.y * cells.z;
  cellList = new int[cellListLen];

  for(int i = natom; i < cellListLen; i ++) cellList[i] = -1;
  for(int i = 0; i < natom; i ++) {
    //tansfor coordination to 0 ~ regInfo.xprd/yprd/zprd
    rs.x = atom_init[i].coord.x - regInfo.xlo;
    rs.y = atom_init[i].coord.y - regInfo.ylo;
    rs.z = atom_init[i].coord.z - regInfo.zlo;
    cc.x = int (rs.x / cellLength);
    cc.y = int (rs.y / cellLength);
    cc.z = int (rs.z / cellLength);
    if(cc.x >= cells.x) cc.x = cells.x - 1;
    if(cc.x <= 0) cc.x = 0;
    if(cc.y >= cells.y) cc.y = cells.y - 1;
    if(cc.y <= 0) cc.y = 0;
    if(cc.z >= cells.z) cc.z = cells.z - 1;
    if(cc.z <= 0) cc.z = 0;
    int c = cc.z * cells.x * cells.y + cc.y * cells.x + cc.x + natom;
    if(c > cellListLen) { printf("BUILDING CELLIST ERROR!\n"); return;}
    cellList[i] = cellList[c];
    cellList[c] = i;
  }
  // to find the nearest atom close to each ctesf[i]
  for(int i = 0; i < ntesf; i ++){
    rs.x = ctesf[i].x - regInfo.xlo;
    rs.y = ctesf[i].y - regInfo.ylo;
    rs.z = ctesf[i].z - regInfo.zlo;
    m1v.x = int (rs.x / cellLength);
    m1v.y = int (rs.y / cellLength);
    m1v.z = int (rs.z / cellLength);
    if(m1v.x >= cells.x) m1v.x = cells.x - 1;
    if(m1v.x <= 0) m1v.x = 0;
    if(m1v.y >= cells.y) m1v.y = cells.y - 1;
    if(m1v.y <= 0) m1v.y = 0;
    if(m1v.z >= cells.z) m1v.z = cells.z - 1;
    if(m1v.z <= 0) m1v.z = 0;

    for(int m1x = -1; m1x <=1; m1x ++){
      for(int m1y = -1 ; m1y <=1; m1y ++){
        for(int m1z = -1; m1z <=1; m1z ++){
          m2v.x = m1v.x + m1x;
          m2v.y = m1v.y + m1y;
          m2v.z = m1v.z + m1z;
          //-------apply periodic boundary condition-----        
          if(m2v.x >= cells.x){ m2v.x = 0; shift.x = regInfo.xprd; } 
          else if(m2v.x < 0){ m2v.x = cells.x - 1; shift.x = -regInfo.xprd; }
          if(m2v.y >= cells.y){ m2v.y = 0; shift.y = regInfo.yprd; }
          else if(m2v.y < 0){ m2v.y = cells.y - 1; shift.y = -regInfo.yprd; }
          if(m2v.z >= cells.z){ m2v.z = 0; shift.z = regInfo.zprd; }
          else if(m2v.z < 0){ m2v.z = cells.z - 1; shift.z = -regInfo.zprd; }
          //--------end of periodic boundary condition-------
          int m2 = m2v.z * cells.x * cells.y + m2v.y * cells.x + m2v.x + natom;
          for(int j2 = cellList[m2]; j2 >= 0; j2 = cellList[j2]) {
            dr.x = atom_init[j2].coord.x - ctesf[i].x;
            dr.y = atom_init[j2].coord.y - ctesf[i].y;
            dr.z = atom_init[j2].coord.z - ctesf[i].z;
            dr.x -= shift.x;
            dr.y -= shift.y;
            dr.z -= shift.z;              
            double rr = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
            if(rr < rCut * rCut){
              ctesf[i] = atom_init[j2].coord;
              Vect norm = normtesf[i].v;
              NormalizeV(norm);
              Vect d = VMul(norm, sqrt(3.0)/6.*lat_a);
              ctesf[i] = VAdd(ctesf[i], d);
              j2 = -1;
              m1x = 2;
              m1y = 2;
              m1z = 2;
            }
          }
        }
      }
    }
  }
  //release memory
  if(nebr != NULL) free(nebr);
  if(cellList != NULL) free(cellList);
}
/*------------------------------------------------------*/
int CreateDisl::IsNumPar(char *par)
{
  int n = strlen(par);
  for(int i = 0; i < n; i ++){
    if(par[i] == '/' && (i == 0 || i == n-1) ){
      return 0;
    }
    else if(par[i] == '-' && i != 0){
      return 0;
    }
    else if(isdigit(par[i]) == 0 && par[i] != '.' && par[i] != '/' && par[i] != '-')
      return 0; //non-numeric
  }
  return 1; //numeric
}
/*------------------------------------------------------*/
void CreateDisl::IndexTransform(Direct &dir)
{
  //transform Miller/Miller-Baravis indexes to a vector in the Cartesian coordinate system
  if(indexType == 1){
    dir.v.x = lat_a * dir.il;
    dir.v.y = lat_a * dir.im;
    dir.v.z = lat_a * dir.in;
  }
  else if(indexType == 2){
    dir.v.x = 1.5 * lat_a * dir.h;
    dir.v.y = sqrt(3.)/2. * lat_a * (dir.h + 2.*dir.k);
    dir.v.z = lat_c * dir.l;
  }
}
/*------------------------------------------------------*/
void CreateDisl::ReadIndex(int indexType, Direct &x, char **items, int nitem, int rowNum)
{
  if(indexType == 1){
    if(nitem != 4){
      printf("The Miller index in line %d is NOT consistent with the parameter list previously set\n", rowNum);
      exit(1);
    } 
    else {
      if(!IsNumPar(items[1]) || !IsNumPar(items[2]) || !IsNumPar(items[3])) {
        printf("non-numeric parameter appears in line %d\n", rowNum);
        exit(1);
      }
      x.il = atof(items[1]);
      x.im = atof(items[2]);
      x.in = atof(items[3]);
      if(x.il == 0 && x.im == 0 && x.in == 0){
        printf("The Miller indexes CANNOT all be set to zeros in line %d\n", rowNum);
        exit(1);
      }
    }
  }
  else if(indexType == 2){
    if(nitem != 5){
      printf("The Miller-Bravais index in line %d is NOT consistent with the parameter list previously set\n", rowNum);
      exit(1);
    } 
    else {
      if(!IsNumPar(items[1]) || !IsNumPar(items[2]) || !IsNumPar(items[3]) || !IsNumPar(items[4])) {
        printf("non-numeric parameter appears in line %d\n", rowNum);
        exit(1);
      }
      x.h = atof(items[1]);
      x.k = atof(items[2]);
      x.i = atof(items[3]);
      x.l = atof(items[4]);
      if(x.h == 0 && x.k == 0 && x.i == 0 && x.l == 0){
        printf("Error: the Miller-Bravais indexes CANNOT all be set to zeros in line %d\n", rowNum);
        exit(1);
      }
      if(x.i != -(x.h + x.k)){
        printf("Error: i != -(h + k) in the Miller-Bravais indexes provided in line %d\n", rowNum);
        exit(1);
      }
    }
  }
  else {
    printf("Error: Wrong index type label for the Miller/Miller-Bravais indexes when reading line %d\n", rowNum);
    exit(1);
  }
}
/*------------------------------------------------------*/
void CreateDisl::ReadLatConst(int indexType, char **items, int nitem, int rowNum)
{
  if(indexType == 1){
    if(nitem != 2){
      printf("Error: the lattice constant setup in line %d is NOT consistent with the cubic lattice\n", rowNum);
      exit(1);
    }
    else{
      if(!IsNumPar(items[1])){
        printf("non-numeric parameter appears in line %d\n", rowNum);
        exit(1);
      }
      lat_a = atof(items[1]);
      lat_c = 0;
    }
  }
  else if(indexType == 2){
    if(nitem != 3){
      printf("Error: the lattice constant setup in line %d is NOT consistent with the HCP lattice\n", rowNum);
      exit(1);
    }
    else{
      if(!IsNumPar(items[1]) || !IsNumPar(items[2])){
        printf("non-numeric parameter appears in line %d\n", rowNum);
        exit(1);
      }
      lat_a = atof(items[1]);
      lat_c = atof(items[2]);
    }
  }
  else {
    printf("Error: the lattice constant setup in line %d is incorrect\n", rowNum);
    exit(1);
  }
}
/*------------------------------------------------------*/
void CreateDisl::InitDirect(Direct &d)
{
  d.il  = d.im  = d.in  = 0;
  d.h   = d.k   = d.i   = d.l  = 0;
  d.v.x = d.v.y = d.v.z = 0;
}
/*------------------------------------------------------*/
/*------------------------------------------------------*/
/*------------------------------------------------------*/
/*------------------------------------------------------*/
/*------------------------------------------------------*/
/*------------------------------------------------------*/
/*------------------------------------------------------*/
/*------------------------------------------------------*/
/*------------------------------------------------------*/
/*------------------------------------------------------*/
