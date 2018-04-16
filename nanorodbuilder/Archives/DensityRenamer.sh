#!/bin/bash
# Script for moving files
shopt -s globstar
	
orig='4dense'
new='2dense'
for f in */ ; do
	echo $f
	mv $f ${f/4dense/$new}
done
for f in * ; do
	echo ${f/4dense/$new}
	if [[ $f = *'4dense'* ]]; then
		mv $f ${f/4dense/$new}
	fi
done
for f in **/* ; do
	echo ${f/4dense/$new}
	if [[ $f = *'4dense'* ]]; then
		mv $f ${f/4dense/$new}
	fi
done

#for f in *.* ; do
#  echo "Calling file $f"
#  name=${f%'.prmtop'}
#  newPath=$gpuPath$name/
#  mkdir $newPath
#  mv $name.* $newPath
#
#done
