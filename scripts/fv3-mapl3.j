#!/usr/bin/env bash

# environment
ARCH=$(uname -s)
FINDPATH="readlink -f"
if [ $ARCH == Darwin ]; then
    FINDPATH=realpath
fi
GEOSBIN=$(dirname $($FINDPATH $0))
GEOSDIR=$(dirname $GEOSBIN)
GEOSLIB=$GEOSDIR/lib
source $GEOSBIN/g5_modules.sh
RUN_CMD="$GEOSBIN/esma_mpirun --n "

# generate config files for un
NX=1
NY=6
IM=48
LM=181
DT=1200
python $GEOSBIN/fv3-mapl3/gen-configs.py $GEOSLIB $NX $NY $IM $LM $DT

# run GEOS.x
$RUN_CMD $((NX*NY)) ./install/bin/GEOS.x cap.yaml |& tee log.run
