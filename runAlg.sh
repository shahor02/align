#!/bin/bash


tarball="algGlo.tar.gz"

if [ -e ${tarball} ]; then
    tar -xvzf ${tarball}
    make
fi

if [ $# -gt 0 ] ; then
    aliroot -b -q $1
fi



