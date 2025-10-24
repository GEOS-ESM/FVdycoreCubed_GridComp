#!/usr/bin/env bash

# environment
ARCH=$(uname -s)
if [ $ARCH == Darwin ]; then
    FINDPATH=realpath
else
    FINDPATH="readlink -f"
fi
RUN_CMD="mpirun -prepend-rank --n "
GEOSBIN=$(dirname $($FINDPATH $0))
GEOSLIB=$GEOSBIN/../lib
source $GEOSBIN/g5_modules.sh

# generate config files for un
NX=1
NY=6
IM=48
LM=181
DT=1200
python $GEOSBIN/fv3-mapl3/gen-config.py $GEOSLIB $NX $NY $IM $LM $DT

# run GEOS.x
$RUN_CMD $((NX*NY)) ./install/bin/GEOS.x cap.yaml |& tee log.run
