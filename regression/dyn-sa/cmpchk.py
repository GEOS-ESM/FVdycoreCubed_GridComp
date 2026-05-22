#!/usr/bin/env python

import numpy as np
import netCDF4 as nc4

BASELINE_VERSION = "v3.0.0-rc.13"
BASELINE_DIR = f"/home/pchakrab/workspace/runs/fv3sa/DynCore/{BASELINE_VERSION}/C48L181/baseline"

def alias(name):
    match name:
        case "U":
            alias = "UV_1"
        case "V":
            alias = "UV_2"
        case _:
            alias = name
    return alias

for state in ["internal"]:
    print("\nSTATE:", state)
    baseline = f"{BASELINE_DIR}/one-hour/fvcore_{state}_checkpoint.nc"
    current = f"checkpoints/last/DYNsa_{state}.nc"
    with nc4.Dataset(baseline) as bas, nc4.Dataset(current) as cur:
        for var in bas.variables:
            if alias(var) in cur.variables:
                bas_var = bas[var][:]
                cur_var = cur[alias(var)][0]
                diff = np.abs(bas_var- cur_var)
                diffnorm = np.linalg.norm(diff)
                print(f" {var:<5} diff(Linf, L2): {np.max(diff)}, {diffnorm}")
