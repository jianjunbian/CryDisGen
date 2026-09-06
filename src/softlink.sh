#!/bin/bash
#---------------------------------------------------
# This script is used to create a soft link to CryDisGen
# in a specified directory
#----------------------------------------------------

if [ $# != 1 ]; then
  echo "USAGE: softlink exe_file"
	exit
fi
if [ ! -e $1 ]; then
  echo "file $1 does not exit"
	exit
fi

dir=`pwd`
file=`basename $1`

src=$dir/$file
loc=~/bin/GenerateDislocation

if [ ! -d $loc ]; then
  echo "The destination directory does not exist!" 
	echo "Please create it first."
	echo " '$loc' "
	exit
fi

des=$loc/$file

ln -s $src $des
echo "soft link is created in $des"
