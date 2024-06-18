# Fortran serialization

We use Serialbox to extract the data around the code we are porting. The methodology goes as follows:

- Annotate files with Serialbox directives, and save them following are XXX.SER with XXX the original filename
- Run `generate.sh` on the file turns them into proper F90 fortran, replacing the original code
- Run simulation, data is dumped into `serialized_data` folder in the scratch directory. The date is dumped as pure binary `.dat` files
- Run `convert_raw_to_netcdf.sh` (see below) to transform the binary to NetCDF (format ingested by the Translate test)

Code lives both in the grid comp and @fvcore on the `dsl/serialbox` branch, which is sync'ed with `dsl/develop`.

# Translate test formating

The translate test system ingest NetCDF formats with a proper naming convention.
 - The savepoint have to be named `XXX-In.nc` (input read) and `XXX-In.nc` (output to compare to), with XXX the name of the translate test (e.g `class TranslateXXX)
 - All buffers in inputs need to be present, _including output buffers_. Inputs is a misnommer, it should read "call signature"
 - All the name of the data (e.g. YYY in the directive `!$ser data YYY=ppp`) are the variable name expected in the Translate `self.data_vars[]` mapping

The simulation will outputs raw .dat files. To format them into proper NetCDF use the `convert_raw_to_netcdf.sh` script by tuning the variable in the script:
 - `EXPERIMENT_SCRATCH_DIR`: experiment scratch directory created by `gcm_run.j`
 - `OUTPUT_NETCDF_DIR`: directory where the NetCDF is going to be dumped
 - `LOG_RANK0_FILE`: `stdout` log dump of rank 0 to a file. Use `mpirun --output-filename logs` to make sure the standard out is written by mpi.

 The script will extract the namelist from the log and turn the raw data to NetCDF using it.
