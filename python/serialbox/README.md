# Translate data generation

We use Serialbox to extract the data around the code we are porting. The methodology goes as follows:

- Annotate files with Serialbox directives, and save them following are XXX.SER with XXX the original filename
- Run `generate.sh` on the file turns them into proper F90 fortran, replacing the original code
- Run simulation, data is dumped into `serialized_data` folder. The date is dumped as pure binary `.dat` files
- Run `ndsl-serialbox_to_netcdf` tool (from NDSL middleware) to transform the binary to NetCDF

Code lives both in the grid comp and @fvcore on the `dsl/serialbox` branch, which is sync'ed with `dsl/develop`.
