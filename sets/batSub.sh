#!/bin/bash

Usage() {
    echo 'Usage: mrgfl.sh [-j <JDL=algJDL50.jdl>] inpFile1 inpFile2 ....'
    exit 1
}

if [ $# -lt 1 ] ; then Usage ;fi

jdl="algJDL50.jdl"

if [ "$1" == "-j" ] ; then
   shift 1 ;
   if [ $# -lt 2 ] ; then Usage ;fi
   jdl=$2
fi

while [ $# -gt 0 ] ; do
    coll=$1
    shift 1 ;
    extension="${coll##*.}"
    if [ "$extension" != "xml" ] ; then
       echo "argument $coll has not .xml extension, skip"
       continue
    fi
    outDir=$(basename "$coll")
    outDir="${outDir%.*}"
    echo 'submitting ' $jdl $coll $outDir
    submit $jdl $coll $outDir
    ps
done

