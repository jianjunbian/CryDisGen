#!/bin/bash

#normal=("0  1  1" "0 -1  1" " 1  0  1" "-1  0  1" "1  1  0" "-1  1  0")
#lx=("1  0  0" "1  0  0" " 0  1  0" " 0  1  0" "0  0  1" " 0  0  1")
#ly=("0  1 -1" "0  1  1" "-1  0  1" "-1  0 -1" "1 -1  0" " 1  1  0")
normal=(" 1  1  1" " 1  1 -1" " 1 -1  1" " 1 -1 -1" "-1 -1 -1" "-1 -1  1" "-1  1 -1" "-1  1  1") 
   lxt=(" 1 -1  0" "-1  1  0" "-1 -1  0" " 1  1  0" " 1 -1  0" "-1  1  0" "-1 -1  0" " 1  1  0")
   lyt=(" 1  1 -2" " 1  1  2" " 1 -1 -2" " 1 -1  2" "-1 -1  2" "-1 -1 -2" "-1  1  2" "-1  1 -2")

echo "" > para_file_array 
echo "Region"     >> para_file_array
echo "-x 1 0 0"   >> para_file_array
echo "-y 0 1 0"   >> para_file_array
echo "-z 0 0 1"   >> para_file_array
echo "-lat 3.615" >> para_file_array
echo "-pbc 1 1 1" >> para_file_array
echo " " >> para_file_array

nx=3   #index of disloation loops along x-axis is in the range: -$nx ~ $nx 
ny=3   #index of disloation loops along y-axis is in the range: -$nx ~ $nx 
nz=3   #index of disloation loops along z-axis is in the range: -$nx ~ $nx 

dx=50  #shift of dislocation center along x-axis
dy=50  #shift of dislocation center along y-axis
dz=50  #shift of dislocation center along z-axis

#ntot=`echo "(2*$nx+1)*(2*$ny+1)*(2*$nz+1);scale=1" | bc`
ntot=`echo "$nx*$ny*$nz;scale=1" | bc`
echo "SFTetrahedron $ntot" >> para_file_array
echo "total loops: $ntot" 

RANDOM=13
for((i=-1;i<=1;i++)); do
  for((j=-1;j<=1;j++)); do
	  for((k=-1;k<=1;k++)); do
		  
		  xc=`echo "$dx*$i;scale=2" | bc`
		  yc=`echo "$dx*$j;scale=2" | bc`
		  zc=`echo "$dx*$k;scale=2" | bc`
      
      nid=$(($RANDOM % 8))
			echo "Tetrahedron           "    >> para_file_array
			echo "-tc     $xc $yc $zc   "    >> para_file_array
			echo "-tr     20            "    >> para_file_array
			echo "-tnorm  ${normal[$nid]} "  >> para_file_array
			echo "-lxt    ${lxt[$nid]}     " >> para_file_array
			echo "-lyt    ${lyt[$nid]}     " >> para_file_array
			echo " "                         >> para_file_array

		done
  done
done



