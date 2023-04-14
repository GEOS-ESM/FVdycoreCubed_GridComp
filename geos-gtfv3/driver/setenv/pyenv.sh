VENV_DIR="${VENV_DIR:-}"
GTFV3_DIR="${GTFV3_DIR:-}"

# Sanitize
if [ -z "$VENV_DIR" ]; then
    echo "Must declare VENV_DIR, directory to build/load virtual env for @gtfv3"
fi
if [ -z "$GTFV3_DIR" ]; then
    echo "Must declare GTFV3_DIR, directory to locate @gtfv3"
fi

# Build or load
if [ -f "$VENV_DIR/bin/activate" ]; then
    echo "-=: Virtual env already built @ $VENV_DIR, sourcing... :=-"
    source $VENV_DIR/bin/activate
    echo "-=: Virtual env sourced. :=-"
else
    echo "-=: Virtual env missing, building @ $VENV_DIR :=-"
    python3 -m venv $VENV_DIR
    source $VENV_DIR/bin/activate
    pip3 install -U pip wheel
    cd $GTFV3_DIR
    git submodule init
    git submodule update
    pip3 install -r requirements_dev.txt
    pip3 install cupy-cuda11x
    pip3 install cffi
    cd -
    echo "-=: Virtual env sourced. :=-"
fi

export PYTHONPATH=$PYTHONPATH:$GTFV3_DIR/../geos-gtfv3/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GTFV3_DIR/../geos-gtfv3/driver/build/