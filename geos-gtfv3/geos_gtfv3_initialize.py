import f90nml
from typing import Tuple
from datetime import datetime

import pace.dsl.stencil
import pace.util as util
from pace.stencils.testing.grid import Grid
from pace.util.grid import DampingCoefficients, GridData, MetricTerms

from fv3core import DynamicalCore
from fv3core._config import DynamicalCoreConfig
from fv3core.testing import TranslateFVDynamics

from geos_gtfv3_debug import write_sum_of_vars


#-GLOBAL-variables-start-
dycore_config = None
driver_object = None
dycore = None
#-GLOBAL-variables-end-


def setup_namelist(npx, npy, npz, ntiles, layout) -> DynamicalCoreConfig:
    '''
    '''
    global dycore_config

    if dycore_config is None:
        namelist = f90nml.read("input.nml")
        dycore_config = DynamicalCoreConfig.from_f90nml(namelist)
        assert dycore_config.npx == npx, '{} =/ {} (from GEOS)'.format(dycore_config.npx, npx)
        assert dycore_config.npy == npy, '{} =/ {} (from GEOS)'.format(dycore_config.npy, npy)
        assert dycore_config.npz == npz, '{} =/ {} (from GEOS)'.format(dycore_config.npz, npz)
        assert dycore_config.ntiles == ntiles, '{} =/ {} (from GEOS)'.format(dycore_config.ntiles, ntiles)
        assert dycore_config.layout == layout, '{} =/ {} (from GEOS)'.format(dycore_config.layout, layout)

    return dycore_config


def setup_dycore(
        backend, comm,
        npx, npy, npz,
        is_, ie, js, je, isd, ied, jsd, jed,
        state_in) -> Tuple[TranslateFVDynamics, DynamicalCore]:
    '''
    '''
    global driver_object, dycore

    if dycore_config is None:
        raise ValueError('dycore_config has not been initialized')

    if dycore is None:

        # Setup grid
        layout = dycore_config.layout
        partitioner = util.CubedSpherePartitioner(util.TilePartitioner(layout))
        communicator = util.CubedSphereCommunicator(comm, partitioner)
        grid = Grid.from_namelist(dycore_config, comm.Get_rank(), backend)

        # Setup stencils
        stencil_config = pace.dsl.stencil.StencilConfig(
            backend=backend,
            rebuild=False,
            validate_args=True,
        )
        stencil_factory = pace.dsl.stencil.StencilFactory(
            config=stencil_config,
            grid_indexing=grid.grid_indexing,
        )

        # Setup structures for grid data computations
        metric_terms = MetricTerms.from_tile_sizing(
            npx=dycore_config.npx,
            npy=dycore_config.npy,
            npz=dycore_config.npz,
            communicator=communicator,
            backend=backend,
        )

        # Create initial  state
        driver_object = TranslateFVDynamics([grid], dycore_config, stencil_factory)
        state = driver_object.state_from_inputs(state_in)

        # Instantiate DynamicalCore
        dycore = DynamicalCore(
            comm=communicator,
            grid_data=GridData.new_from_metric_terms(metric_terms),
            stencil_factory=stencil_factory,
            damping_coefficients=DampingCoefficients.new_from_metric_terms(metric_terms),
            config=dycore_config,
            phis=state.phis,
        )

    return driver_object, dycore
