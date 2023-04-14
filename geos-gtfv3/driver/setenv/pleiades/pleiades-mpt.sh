SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

module load comp-gcc/9.3.0-TOSS3
module load mpi-hpe/mpt.2.25 mpi-mpt/2.25/gcc/9.3.0-TOSS3
module load cuda/11.0
module use -a $NOBACKUP/sw/modulefiles
module load pchakrab/python/3.8/gcc9.3.0-cuda11.0-mpt2.25-toss3

export MPI_USE_CUDA=true
export MPI_UNBUFFERED_STDIO=true
export PYTHONPATH=$PYTHONPATH:${SCRIPT_DIR}/../..:$NOBACKUP/code/ext/vulcan/fv3core