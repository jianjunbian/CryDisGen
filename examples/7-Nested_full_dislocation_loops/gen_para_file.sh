#!/bin/bash
 slipv=(" 1  1  1" " 1  1  1")
normal=(" 1  1  0" " 1 -1  0") #( " 1  0  1" "-1  0  1" "1  1  0" "-1  1  0")
    lx=("-1  1  0" " 1  1  0") #( " 0  1  0" " 0  1  0" "0  0  1" " 0  0  1")
    ly=(" 0  0  1" " 0  0  1") #( "-1  0  1" "-1  0 -1" "1 -1  0" " 1  1  0")

echo "" > para_file_array 
echo "Region"     >> para_file_array
echo "-x 1 0 0"   >> para_file_array
echo "-y 0 1 0"   >> para_file_array
echo "-z 0 0 1"   >> para_file_array
echo "-lat 3.303" >> para_file_array
echo "-pbc 1 1 1" >> para_file_array
echo " " >> para_file_array


nx=0   #index of disloation loops along x-axis is in the range: -$nx ~ $nx 
ny=0   #index of disloation loops along x-axis is in the range: -$nx ~ $nx 
nz=2   #index of disloation loops along x-axis is in the range: -$nx ~ $nx 

dx=0   #shift of dislocation center
dy=0
dz=30

ntot=`echo "(2*$nx+1)*(2*$ny+1)*(2*$nz+1);scale=1" | bc`
echo "Loop $ntot" >> para_file_array
echo "total loops: $ntot" 

n=0
for((i=0;i<=0;i++)); do
  for((j=0;j<=0;j++)); do
	  for((k=-$nz;k<=$nz;k++)); do
		  
		  xc=`echo "$dx*$i;scale=2" | bc`
		  yc=`echo "$dy*$j;scale=2" | bc`
		  zc=`echo "$dz*$k;scale=2" | bc`
			echo $zc

			n=$((n+1))
			nid=$((n % 2))

			echo "Dislocation           "   >> para_file_array
			echo "-bv     ${slipv[$nid]}"   >> para_file_array
			echo "-bs     1/2           "   >> para_file_array
			echo "-dnorm  ${normal[$nid]} " >> para_file_array
			echo "-dc     $xc $yc $zc   "   >> para_file_array
			echo "-dr     25            "   >> para_file_array
			echo "-lx     ${lx[$nid]}     " >> para_file_array
			echo "-ly     ${ly[$nid]}     " >> para_file_array
			echo "-ndisl  8             "   >> para_file_array
			echo " "                        >> para_file_array

		done
  done
done



