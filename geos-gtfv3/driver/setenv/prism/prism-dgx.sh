SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

module use -a /att/nobackup/pchakrab/sw/modulefiles
module load gcc/9.2.0
module load pchakrab/cuda/11.5
module load pchakrab/ucx/1.11.2/gcc9.2.0-cuda11.5
module load pchakrab/mpi/openmpi/4.0.6/gcc9.2.0-cuda11.5
module load pchakrab/python/3.8/gcc9.2.0-cuda11.5-ompi4.0.6

export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$NOBACKUP/sw/boost_1_78_0
export PYTHONPATH=${PYTHONPATH}:${SCRIPT_DIR}/../..:$NOBACKUP/code/ext/vulcan/fv3core
