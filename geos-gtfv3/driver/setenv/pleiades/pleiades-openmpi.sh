SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
module load gcc/9.3
module load cuda/11.0
module use -a $NOBACKUP/sw/modulefiles
module load pchakrab/mpi/openmpi/4.0.6/gcc9.3.0-cuda11.0
module load pchakrab/python/3.8/gcc9.3.0-cuda11.0-openmpi4.0.6-toss3

export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$NOBACKUP/sw/boost_1_78_0
export PYTHONPATH=$PYTHONPATH:${SCRIPT_DIR}/../..:$NOBACKUP/code/ext/vulcan/fv3core
