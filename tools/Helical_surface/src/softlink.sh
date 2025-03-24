#!/bin/bash
if [ $# != 1 ]; then
  echo "USAGE: softlink exe_file"
	exit
fi
if [ ! -e $1 ]; then
  echo "file $1 does not exit"
	exit
fi

file=`basename $1`

src_dir=`pwd`
src=$src_dir/$file

des_dir=~/bin/GenerateDislocation
if [ ! -d $des_dir ]; then
  echo "target directory does not exist"
	exit
fi
des=$des_dir/$file

ln -s $src $des
echo "soft link is created in $des"
