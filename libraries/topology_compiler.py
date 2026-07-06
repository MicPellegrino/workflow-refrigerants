import os
import subprocess
import re
import numpy as np

# Ugly but necessary...
EPS_MACHINE = np.finfo(float).eps

# Oldest admissible GMX version
OLDEST_GMX_VER = 2024

# No. lines to keep when running 'insert-molecules'
END_LINES_BUFFER = 20

# Conversion kg/m^3 -> u/nm^3
DENSITY_SI_2_GMX = 0.602214


def test_gromacs_availability(gmx_bin='gmx',
                              oldest_gmx_ver=OLDEST_GMX_VER):

    error_string_bin_not_found = (
    "ERROR: GROMACS binaries not found, or returned an error!\n"
    "Make sure to run 'source <path-to-gromacs-bin>/GMRC'\n"
    "or specify the full path to the binary as input."
    )
    try:
        subprocess.run([gmx_bin], check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(error_string_bin_not_found)

    std_temp = subprocess.run([gmx_bin, "-version"], capture_output=True, text=True)
    version_string = std_temp.stdout + std_temp.stderr

    match = re.search(r'GROMACS - gmx, (\d{4})', version_string)
    if match:
        gmx_version = int(match.group(1))
        if gmx_version < oldest_gmx_ver :
            warning_string_old_version = (
            f"WARNING: Old GROMACS version detected ({gmx_version})!\n"
            f"We recommend using GROMACS {oldest_gmx_ver} or later."
            )
            print(warning_string_old_version)
    else:
        print("WARNING: Could not determine GROMACS version")

test_gromacs_availability.__doc__=f"""
    Input:
    - gmx_bin ('gmx'): Path to GROMACS binary file;
    - oldest_gmx_ver ({OLDEST_GMX_VER}): Oldest GROMACS version compatible with the library.
    """


def run_x2top(gro_file, 
              ff_folder, 
              name=None, 
              top_file=None, 
              flags="", 
              gmx_bin='gmx'):

    assert gro_file.split(".")[-1]=='gro', "Please provide a '.gro' file as input configuration!"

    ff_folder_stem = ff_folder[:-3] if ff_folder.endswith('.ff') else ff_folder

    # BUG!
    # Test the case where the input file is in "../<folder>"!
    # This affect all other calls to .join(...)

    if name == None :
        name = ''.join(gro_file.split(".")[:-1])
    if top_file == None :
        top_file = ''.join(gro_file.split(".")[:-1])+".top"
    
    cmd = f"{gmx_bin} x2top -f {gro_file} -ff {ff_folder_stem} -name {name} -o {top_file} {flags}"
    os.system(cmd)

run_x2top.__doc__=f"""
    Input:
    - gro_file: GROMACS configuration file (.gro) with molecular coordinates;
    - ff_folder: Folder containing interatomic potential definitions;
    - name (None): Name of the molecule in the output .top file (default name: conf. file without extension);
    - top_file (None): Name of the output .top file (default name: conf. file with .top extension);
    - flags (""): Any additional flag to pass to 'gmx x2top';
    - gmx_bin ('gmx'): Path to GROMACS binary file.
    """


def select_top_lines(top_file):

    with open(top_file, 'r') as f:
        top_lines = f.readlines()

    top_lines = [line for line in top_lines if not (line.lstrip().startswith(';') or line.lstrip().startswith('#include'))]

    skip_sections = {'system', 'molecules'}
    itp_lines = []
    skip = False
    for line in top_lines:
        stripped = line.strip()
        if stripped.startswith('[') and stripped.endswith(']'):
            section = stripped[1:-1].strip().lower()
            skip = section in skip_sections
        if not skip:
            itp_lines.append(line)

    return itp_lines


def read_charges(charge_list_file):

    charges_dict = dict()

    tot_charge = 0

    with open(charge_list_file, 'r') as f:
        for line in f:
            if line.lstrip().startswith(';') or not line.strip():
                continue
            parts = line.split()
            atom_name, charge = parts[1], float(parts[2])
            tot_charge += charge
            if atom_name in charges_dict:
                raise ValueError(f"Duplicate atom name '{atom_name}' found in charges file")
            charges_dict[atom_name] = charge

    if np.abs(tot_charge)>EPS_MACHINE :
        print(f"WARNING: Total charge does not add to zero (q_tot = {tot_charge} e).")

    return charges_dict


def fix_topology_file(topology_file,
                      charge_list_file):

    charges_dict = read_charges(charge_list_file)

    with open(topology_file, 'r') as f:
        lines = f.readlines()

    result = []
    in_atoms = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('[') and stripped.endswith(']'):
            in_atoms = stripped == '[ atoms ]'
        elif in_atoms and not stripped.startswith(';') and stripped:
            parts = line.split()
            if parts[4] in charges_dict:
                new_charge = f"{charges_dict[parts[4]]:.5f}"
                line = line[:line.index(parts[6])] + new_charge + line[line.index(parts[6]) + len(parts[6]):]
        result.append(line)

    extension = topology_file.split(".")[-1]
    topology_file_fixed = ''.join(topology_file.split(".")[:-1])+"-fixed."+extension

    with open(topology_file_fixed, 'w') as f:
        f.writelines(result)

    return topology_file_fixed


def create_itp(top_file,
               charge_list_file=None,
               itp_file=None):

    assert top_file.split(".")[-1]=='top', "ERROR: Please provide a '.top' file as input topology!"

    if charge_list_file==None :
        itp_lines = select_top_lines(top_file)
    else :
        top_file_fixed = fix_topology_file(top_file,charge_list_file)
        itp_lines = select_top_lines(top_file_fixed)

    if itp_file == None :
        itp_file = ''.join(top_file.split(".")[:-1])+".itp"

    with open(itp_file, 'w') as f:
        f.writelines(itp_lines)

create_itp.__doc__ = f"""
    Input:
    - top_file: Topology file (.top) of the molecule;
    - charge_list_file (None): List of partial charges for each atom in the molecule (.txt);
    - itp_file: Output .itp file (default name: .top file with .itp extension).
    """


def run_insert_molecules(gro_file_f,
                         gro_file_ci,
                         nmol, 
                         flags="",
                         gro_file_out=None,
                         gmx_bin='gmx'):

    assert gro_file_f.split(".")[-1]=='gro', "ERROR: Please provide a '.gro' file as input configuration!"
    assert gro_file_ci.split(".")[-1]=='gro', "ERROR: Please provide a '.gro' file as input solvant!"

    if gro_file_out==None :
        gro_file_out = ''.join(gro_file_f.split(".")[:-1])+"-"+''.join(gro_file_ci.split(".")[:-1])+".gro"

    cmd = f"{gmx_bin} insert-molecules -f {gro_file_f} -ci {gro_file_ci} -nmol {nmol} -o {gro_file_out} {flags}"
    stderr_lines = []
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
        for line in proc.stderr:
            if line != '\n':
                print(line, end="")
                stderr_lines.append(line)
    tail = stderr_lines[-END_LINES_BUFFER:]

    n_added_mol = None
    for line in tail:
        m = re.search(r"Added (\d+) molecules \(out of (\d+) requested\)", line)
        if m:
            n_added_mol = int(m.group(1))
            break
    if n_added_mol is None:
        raise RuntimeError("Could not parse molecule count from GROMACS output.")

    return n_added_mol

run_insert_molecules.__doc__ = f"""
    Input:
    - gro_file_f: GROMACS .gro file for the substrate/surface;
    - gro_file_ci: GROMACS .gro file for the solvant;
    - nmol: Number of solvant molecules to insert;
    - flags (""): Any additional flag to pass to "gmx insert-molecules";
    - gro_file_out (None): Output configuration file (default name: combine solvant and substrate with .gro extension)
    - gmx_bin ('gmx'): Path to GROMACS binary file.
    Output:
    - n_added_mol: number of solvane molecules added to the configuration.
    """


def extract_moleculetype_name(filepath):

    with open(filepath, 'r') as f:
        content = f.read()
    match = re.search(
        r'\[\s*moleculetype\s*\]\s*(?:;[^\n]*\n\s*)*([^\s;][^\n]*)',
        content
    )
    if not match:
        raise ValueError("No [ moleculetype ] section found.")
    name = match.group(1).split()[0]
    return name


def compile_topology(n_sol,
                     n_sub,
                     ff_itp,
                     sol_itp,
                     sub_snippet,
                     output_top="topology-biphase.top",
                     system_name="Biphase"):

    mol_name_sol = extract_moleculetype_name(sol_itp)
    mol_name_sub = extract_moleculetype_name(sub_snippet)

    sol_itp_path = os.path.realpath(os.path.abspath(os.path.expanduser(sol_itp)))
    ff_itp_path = os.path.realpath(os.path.abspath(os.path.expanduser(ff_itp)))

    with open(output_top, 'w') as fo:
        fo.write(f"""#include "{ff_itp_path}" \n""")
        fo.write( """\n""")
        with open(sub_snippet, 'r') as fi:
            sub_snippet_lines = fi.read()
        fo.write(sub_snippet_lines)
        fo.write( """\n""")
        fo.write(f"""#include "{sol_itp_path}" \n""")
        fo.write( """\n""")
        fo.write( """[ system ]\n""")
        fo.write(f"""{system_name}\n""")
        fo.write( """\n""")
        fo.write( """[ molecules ]\n""")
        fo.write(f"""{mol_name_sub} {n_sub}\n""")
        fo.write(f"""{mol_name_sol} {n_sol}\n""")

compile_topology.__doc__ = f"""
    Input:
    - n_sol: No. of solvant molecules;
    - n_sub: No. of substrate/surface molecules;
    - ff_itp: Path to ff file (header with [ defaults ] definition);
    - sol_itp: Topology of solvant molecule;
    - sub_snippet: Header of substrate topology (e.g. [ atomtypes ], ..., [ position_restraints ], ...);
    - output_top ("topology-biphase.top"): ...;
    - system_name ("Biphase"): ....
    """


if __name__=="__main__":

    test_gromacs_availability()
    run_x2top("HFO-1234zeE.gro", "refrigerants.ff")
    create_itp("HFO-1234zeE.top", charge_list_file="charges.txt")
    n_added_mol = run_insert_molecules("zirconia.gro", "HFO-1234zeE.gro", flags="-try 10 -scale 0.65", nmol=1800)
    compile_topology(n_added_mol, 1600, "refrigerants.ff/refrig_nb.itp", "HFO-1234zeE.itp", "zirconia.txt")