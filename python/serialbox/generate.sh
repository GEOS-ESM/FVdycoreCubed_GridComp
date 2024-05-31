#!/bin/bash

ml SMTStack/2024.04.00

FV_COMP=../../

FILE=$FV_COMP/FV_StateMod.F90
python3 $SERIALBOX_ROOT/python/pp_ser/pp_ser.py $FILE.SER > $FILE
tail -n +2 $FILE > $FILE.swp && mv $FILE.swp $FILE
