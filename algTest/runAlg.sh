#!/bin/bash

export LD_PRELOAD=$LD_PRELOAD:/usr/lib/x86_64-linux-gnu/libjemalloc.so

tarball="algGlo.tar.gz"

if [ -e ${tarball} ]; then
    tar -xvzf ${tarball}
    make
fi

prc='runGloAlgTask.C'
if [ $# -gt 0 ] ; then
#  ddd  aliroot  $1
    prc=$1
fi

aliroot -b -q $prc


