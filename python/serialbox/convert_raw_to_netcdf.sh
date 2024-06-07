EXPERIMENT_SCRATCH_DIR=/home/fgdeconi/work/git/fp/experiment/TBC_C24_L72/scratch
OUTPUT_NETCDF_DIR=/home/fgdeconi/work/git/fp/experiment/TBC_C24_L72/SavepointNetcdf

cp $EXPERIMENT_SCRATCH_DIR/input.nml $EXPERIMENT_SCRATCH_DIR
ndsl-serialbox_to_netcdf $EXPERIMENT_SCRATCH_DIR $OUTPUT_NETCDF_DIR
