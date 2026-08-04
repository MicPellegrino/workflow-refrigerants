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