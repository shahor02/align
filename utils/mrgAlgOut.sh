#!/bin/bash

Usage() {
    echo 'Usage: mrgAlgOut.sh <source_dir1> <source_dir2> ...'
    echo 'Merged output will be stored in source directory'
    exit 1
}

control=mpControlRes.root
stat=mpStatOut.root
mpdat=mpData.root

if [ $# -lt 1 ] ; then Usage ;fi

for dir in "$@"; 
do
#
 dir=`echo $dir | sed "s/\/$//"`
 echo doing $dir
# mrgfl.sh -f ${dir}/$control ${dir}/*/$control
# mrgfl.sh -f ${dir}/$stat ${dir}/*/$stat
 mrgfl.sh -f ${dir}/$mpdat ${dir}/*/$mpdat

#
done

