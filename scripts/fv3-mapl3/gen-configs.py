#!/usr/bin/env python

import os
import sys
import yaml
import jinja2

script_dir = os.path.dirname(__file__)

def gen_cap_yaml(lib_dir, nx, ny):
    with open(os.path.join(script_dir, "cap.j2")) as fin:
        cap_template = jinja2.Template(fin.read())
    cap_yaml = cap_template.render(lib_dir=lib_dir, petcount=nx*ny)
    with open("cap.yaml", "w") as fout:
        fout.write(cap_yaml)

def gen_fv3sa_yaml(nx, ny, im, num_levels, dt):
    with open(os.path.join(script_dir, "fv3sa.j2")) as fin:
        fv3sa_template = jinja2.Template(fin.read())
    fv3sa_yaml = fv3sa_template.render(
        nx=nx, ny=ny,
        im=im, jm=im*6, num_levels=num_levels,
        run_dt=dt)
    with open("fv3sa.yaml", "w") as fout:
        fout.write(fv3sa_yaml)

def gen_input_nml(im, num_levels):
    with open(os.path.join(script_dir, "input.j2")) as fin:
        input_template = jinja2.Template(fin.read())
    input_nml = input_template.render(npx=im+1, npy=im+1, npz=num_levels)
    with open("input.nml", "w") as fout:
        fout.write(input_nml)

def gen_cap_restart_yaml():
    with open("cap_restart.yaml", "w") as fout:
        fout.write('currTime: "1891-03-02T00:00:00"\n')

if __name__ == "__main__":
    lib_dir = sys.argv[1]
    nx = int(sys.argv[2])
    ny = int(sys.argv[3])
    im = int(sys.argv[4])
    num_levels = int(sys.argv[5])
    dt = int(sys.argv[6])

    gen_cap_yaml(lib_dir, nx, ny)
    gen_fv3sa_yaml(nx, ny, im, num_levels, dt)
    gen_input_nml(im, num_levels)
    gen_cap_restart_yaml()
