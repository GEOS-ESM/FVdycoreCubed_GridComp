#!/usr/bin/env python

import numpy as np
import netCDF4 as nc4

BASELINE_VERSION = "GCMv12-rc27p1+Held-Suarez"
BASELINE_DIR = f"/home/pchakrab/workspace/runs/fv3sa/AdvCore/{BASELINE_VERSION}/C12L91/baseline/const-tracer"

def alias(name):
    match name:
        case "U":
            alias = "UV_1"
        case "V":
            alias = "UV_2"
        case _:
            alias = name
    return alias

print("\nDYN")
states = ["internal"]
for state in states:
    print(f"\n{state}")
    baseline = f"{BASELINE_DIR}/one-hour/DYN_{state}_checkpoint.nc"
    current = f"checkpoints/last/DYN_{state}.nc"
    with nc4.Dataset(baseline) as bas, nc4.Dataset(current) as cur:
        for var in bas.variables:
            if alias(var) in cur.variables:
                bas_var = bas[var][:]
                cur_var = cur[alias(var)][0]
                diff = np.abs(bas_var- cur_var)
                diffnorm = np.linalg.norm(diff)
                print(f" {var:<5} diff(Linf, L2): {np.max(diff)}, {diffnorm}")

print("\nADV")
states = ["export"]
for state in states:
    print(f"\n{state}")
    baseline = f"{BASELINE_DIR}/one-hour/ADV_{state}_checkpoint.nc"
    current = f"checkpoints/last/ADV_{state}.nc"
    with nc4.Dataset(baseline) as bas, nc4.Dataset(current) as cur:
        for var in bas.variables:
            if var in cur.variables:
                bas_var = bas[var][:]
                cur_var = cur[var][0]
                diff = np.abs(bas_var- cur_var)
                diffnorm = np.linalg.norm(diff)
                print(f" {var:<5} diff(Linf, L2): {np.max(diff)}, {diffnorm}")
