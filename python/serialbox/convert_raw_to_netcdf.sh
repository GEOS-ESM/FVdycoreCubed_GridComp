EXPERIMENT_SCRATCH_DIR=/home/mad/work/fp/experiments/TBC_C24_L72/scratch
OUTPUT_NETCDF_DIR=/home/mad/work/fp/savepoints

cp $EXPERIMENT_SCRATCH_DIR/input.nml $EXPERIMENT_SCRATCH_DIR/serialized_data
ndsl-serialbox_to_netcdf $EXPERIMENT_SCRATCH_DIR/serialized_data $OUTPUT_NETCDF_DIR
