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
FUN_AREA = lambda theta_0, R_0 : R_0*R_0*(theta_0/(np.sin(theta_0)**2)-1.0/np.tan(theta_0))
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


# TODO: check if this works on a generic .gro file
def parse_atom_line(line):

    resname = line[5:10].strip()
    x = float(line[20:28])
    y = float(line[28:36])
    z = float(line[36:44])
    
    return x, y, z, resname


def carve_gro(input_file,
              atoms_per_molecule,
              carve_condition,
              mol_type='UNK',
              output_file=None):

    # TODO: assertions and deal with case output_file==None

    with open(input_file, "r") as f:
        lines = f.readlines()
 
    title = lines[0]
    n_atoms = int(lines[1].strip())
    atom_lines = lines[2:2 + n_atoms]
    box_line = lines[2 + n_atoms:]
 
    kept_lines = []
 
    i = 0
    while i<n_atoms :

        # for i in range(0, n_atoms, atoms_per_molecule):
        x, y, z, resname = parse_atom_line(atom_lines[i])

        if resname != mol_type :
            kept_lines.append(atom_lines[i])
            i = i+1

        else :
            molecule_lines = atom_lines[i:i + atoms_per_molecule]
 
            passes = True
            for line in molecule_lines:
                x, y, z, resname = parse_atom_line(line)
                if not carve_condition(x, y, z):
                    passes = False
                    break
 
            if passes:
                kept_lines.extend(molecule_lines)

            i += atoms_per_molecule

    print("Total number of new atoms:",len(kept_lines))
 
    with open(output_file, "w") as f:
        f.write(title)
        f.write(f"{len(kept_lines)}\n")
        f.writelines(kept_lines)
        f.writelines(box_line)


if __name__=="__main__":

    lx_0 = 5.26800
    nbox_x = 10
    R0 = 0.5*(nbox_x*lx_0)-5.0
    theta_0 = np.deg2rad(10)
    r0 = FUN_RADIUS(theta_0,R0)
    print("r0 =",r0,"nm")

    extend_substrate("../example-droplet/zirconia.gro",nbox_x,output_file="../example-droplet/zirconia-ext.gro")
    solvate_empty_box("../example-droplet/box-HFO-1234zeE.gro", "../example-droplet/empty.gro", "../example-droplet/solvated.gro")

    cx = 26.34
    cz = 7.27545
    carve_condition = lambda x, y, z : ((x-cx)*(x-cx)+(z-cz)*(z-cz))<=(r0*r0)

    n_atom_per_mol = 9
    carve_gro("../example-droplet/solvated.gro",n_atom_per_mol,carve_condition,output_file="../example-droplet/carved.gro")

    # TODO: Insert molecules to fill the vapour phase

    zupp = 12.5
    zlow = 1.5
    carve_condition = lambda x, y, z : (z>zlow)*(z<zupp)
    carve_gro("../example-droplet/md-droplet/nvt.gro",n_atom_per_mol,carve_condition,output_file="../example-droplet/carved-temp.gro")