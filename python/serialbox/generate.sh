#!/bin/bash

ml SMTStack/2024.04.00

FV_COMP=/home/mad/work/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp

# find ../../ -iname *F90.SER

PPSER_VERB="$SERIALBOX_ROOT/python/pp_ser/pp_ser.py --verbose --ignore-identical -m utils_ppser_kbuff"


FILE=$FV_COMP/FV_StateMod.F90
python3 $PPSER_VERB $FILE.SER > $FILE
tail -n +2 $FILE > $FILE.swp && mv $FILE.swp $FILE

FILE=$FV_COMP/@fvdycore/model/dyn_core.F90
python3 $PPSER_VERB $FILE.SER > $FILE
tail -n +2 $FILE > $FILE.swp && mv $FILE.swp $FILE

FILE=$FV_COMP/@fvdycore/model/fv_dynamics.F90
python3 $PPSER_VERB $FILE.SER > $FILE
tail -n +2 $FILE > $FILE.swp && mv $FILE.swp $FILE

FILE=$FV_COMP/@fvdycore/model/sw_core.F90
python3 $PPSER_VERB $FILE.SER > $FILE
tail -n +2 $FILE > $FILE.swp && mv $FILE.swp $FILE

FILE=$FV_COMP/@fvdycore/model/nh_utils.F90
python3 $PPSER_VERB $FILE.SER > $FILE
tail -n +2 $FILE > $FILE.swp && mv $FILE.swp $FILE
