SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# PACE directories
PACEDIR=/discover/nobackup/pchakrab/code/ext/ai2/pace/6d5bb26
export PYTHONPATH=${PYTHONPATH}:${PACEDIR}/dsl
export PYTHONPATH=${PYTHONPATH}:${PACEDIR}/util
export PYTHONPATH=${PYTHONPATH}:${PACEDIR}/driver
export PYTHONPATH=${PYTHONPATH}:${PACEDIR}/physics
export PYTHONPATH=${PYTHONPATH}:${PACEDIR}/fv3core
export PYTHONPATH=${PYTHONPATH}:${PACEDIR}/stencils

# Directory containing geos-gtfv3.py
export PYTHONPATH=${PYTHONPATH}:${SCRIPT_DIR}/../../..
