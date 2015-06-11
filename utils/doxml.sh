#!/bin/bash
# create in batch XML sets from run list

Usage() {
    echo 'Usage: doxml.sh -d <dir> -p <pass> [-c prefix=coll] run1 [run2 ...]'
    echo 'will crete XML sets for provided run list, as'
    echo 'looking in <dir>/000<run?>/<pass> ...'
    echo 'output to <prefix><run1>.xml <prefix><run2>.xml ...'
    echo 'Example:'
    echo './doxml.sh -d  /alice/data/2015/LHC15a -p cosmics_pass2  -c LHC15a_pass2_ <runs>'
    exit 1
}



if [ $# -lt 1 ] ; then Usage ;fi

dir=''
if [ "$1" == "-d" ] ; then
    shift 1 ;
    if [ $# -lt 2 ] ; then Usage ;fi
    dir="$1"
    shift 1 ;
else
    Usage
fi

pass=''
if [ "$1" == "-p" ] ; then
    shift 1 ;
    if [ $# -lt 2 ] ; then Usage ;fi
    pass="$1"
    shift 1 ;
else
    Usage
fi

if [ $# -lt 1 ] ; then Usage ;fi

pref='coll'
if [ "$1" == "-c" ] ; then
    shift 1 ;
    if [ $# -lt 2 ] ; then Usage ;fi
    pref="$1"
    shift 1 ;
else
    Usage
fi

esdname='AliESDs.root'
if [[ $dir == *cpass* ]] ; then esdname="AliESDs_Barrel.root" ; fi

if [ $# -lt 1 ] ; then Usage ;fi

while [ $# -gt 0 ] ; do
    run=$1; 
    shift 1 ;
    pth=${dir}/000${run}/${pass}
#
    collname=${pref}${run}.xml
    echo doing run $run from path $dir to $collname  
    echo aliroot -b -q CreateDataSet.C\(\"${pth}\",-1,\"${esdname}\",\"${collname}\"\)
    aliroot -b -q CreateDataSet.C\(\"${pth}\",-1,\"${esdname}\",\"${collname}\"\)
#
done
