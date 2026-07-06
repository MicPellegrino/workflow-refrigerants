import os
import subprocess
import re
import numpy as np

# TODO: these variables should be defined in a separate library

# Ugly but necessary...
EPS_MACHINE = np.finfo(float).eps

# Oldest admissible GMX version
OLDEST_GMX_VER = 2024

# No. lines to keep when running 'insert-molecules'
END_LINES_BUFFER = 20

# Conversion kg/m^3 -> u/nm^3
DENSITY_SI_2_GMX = 0.602214

# NB: theta_0 is in radians!
FUN_AREA = lambda theta_0, R_0 : R_0*R_0*(theta_0/np.sin(theta_0)**2-1.0/np.tan(theta_0))
FUN_RADIUS = lambda theta_0, R_0 : np.sqrt(FUN_AREA(theta_0,R_0)/np.pi)

# TODO: test_gromacs_availability() in topolgy_compiler.py should be a global function

def extend_substrate(input_file,
                     nx,
                     ny=1,
                     nz=1,
                     output_file=None,
                     flags="",
                     gmx_bin='gmx'):

    assert input_file.split(".")[-1]=='gro', "ERROR: Please provide a '.gro' file as input configuration!"

    # BUG
    if output_file == None :
        output_file = ''.join(input_file.split(".")[:-1])+"-ext.gro"

    cmd = f"{gmx_bin} genconf -nbox {nx} {ny} {nz} -f {input_file} -o {output_file} {flags}"
    os.system(cmd)


# TODO: the empty box should be created by the function, given its size
def solvate_empty_box(input_solvant,
                      input_empty,
                      output_file="solvated.gro",
                      flags="",
                      gmx_bin='gmx'):

    # TODO: assertions...

    cmd = f"{gmx_bin} solvate -cp {input_empty} -cs {input_solvant} -o {output_file} {flags}"
    os.system(cmd)


if __name__=="__main__":

    # extend_substrate("../example/zirconia.gro",5,output_file="../example/zirconia-ext.gro")
    solvate_empty_box("../example/box-HFO-1234zeE.gro", "../example/empty.gro", "../example/solvated.gro")

    theta_0 = np.deg2rad(10)
    R0 = 10.0
    print("r0 =", FUN_RADIUS(theta_0,R0),"nm")
