#!/bin/bash

normal=("1  1  1" "1  1 -1" "-1  1  1" "-1  1 -1")
lx=(" 1 -1  0" " 1 -1  0" " 1  1  0" " 1  1  0")
ly=(" 1  1 -2" "-1 -1 -2" "-1  1 -2" " 1 -1 -2")

echo "" > para_file_array 
echo "Region"     >> para_file_array
echo "-x 1 0 0"   >> para_file_array
echo "-y 0 1 0"   >> para_file_array
echo "-z 0 0 1"   >> para_file_array
echo "-lat 3.301" >> para_file_array
echo "-pbc 1 1 1" >> para_file_array
echo " " >> para_file_array


nx=1   #index of disloation loops along x-axis is in the range: -$nx ~ $nx 
ny=1   #index of disloation loops along x-axis is in the range: -$nx ~ $nx 
nz=1   #index of disloation loops along x-axis is in the range: -$nx ~ $nx 

dx=60  #shift of dislocation center
dy=60
dz=60

ntot=`echo "(2*$nx+1)*(2*$ny+1)*(2*$nz+1);scale=1" | bc`
echo "Loop $ntot" >> para_file_array
echo "total loops: $ntot" 

RANDOM=13

for((i=-$nx;i<=$nx;i++)); do
  for((j=-$ny;j<=$ny;j++)); do
	  for((k=-$nz;k<=$nz;k++)); do
		  
		  xc=`echo "$dx*$i;scale=2" | bc`
		  yc=`echo "$dx*$j;scale=2" | bc`
		  zc=`echo "$dx*$k;scale=2" | bc`

			nid=$(($RANDOM % 4))

			echo "Dislocation           "   >> para_file_array
			echo "-bv     ${normal[$nid]} " >> para_file_array
			echo "-bs     1/2           "   >> para_file_array
			echo "-dnorm  ${normal[$nid]} " >> para_file_array
			echo "-dc     $xc $yc $zc   "   >> para_file_array
			echo "-dr     20            "   >> para_file_array
			echo "-lx     ${lx[$nid]}     " >> para_file_array
			echo "-ly     ${ly[$nid]}     " >> para_file_array
			echo "-ndisl  6             "   >> para_file_array
			echo " "                        >> para_file_array

		done
  done
done



