#!/bin/bash

for folder in "$@"
do
	group=`basename $folder`
	for file in $folder/*.hsq
	do
		fn=${file##*/}
		pheno=${fn%.*}
		est=$(awk '$1=="V(G)/Vp"{print $2,$3}' $file)
		if [ -n "$est" ]; then
			echo $group $pheno $est
		fi
	done
done
