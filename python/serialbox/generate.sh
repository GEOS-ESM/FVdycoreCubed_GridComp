#!/bin/bash

ml SMTStack/2024.04.00

FV_COMP=/home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp

# find ../../ -iname *F90.SER

FILE=$FV_COMP/FV_StateMod.F90
python3 $SERIALBOX_ROOT/python/pp_ser/pp_ser.py $FILE.SER > $FILE
tail -n +2 $FILE > $FILE.swp && mv $FILE.swp $FILE

FILE=$FV_COMP/@fvdycore/model/dyn_core.F90
python3 $SERIALBOX_ROOT/python/pp_ser/pp_ser.py $FILE.SER > $FILE
tail -n +2 $FILE > $FILE.swp && mv $FILE.swp $FILE
