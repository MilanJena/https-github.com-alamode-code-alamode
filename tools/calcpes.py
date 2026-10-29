import argparse

import numpy as np
from alamode_constants import BOHR_IN_ANGSTROM as Bohr
from alamode_constants import RYD_IN_EV as Ry
from ase import Atoms
from ase.io import read

from fcsio.reader import ForceConstantParser
from taylor import TaylorExpansionPotential


def get_original_structure(fname_in, format):
    """
    Read the original structure from a file.

    Parameters:
    - fname_in (str): The file name of the structure file.
    - format (str): The format of the structure file. If None, the format is auto-detected.

    Returns:
    - structure: The structure object read from the file.
    """
    if format is None:
        structure = read(fname_in)
    else:
        structure = read(fname_in, format=format)

    return structure


def print_energies(energies, unit="rydberg"):
    """
    Print the calculated potential energies.

    Parameters:
    - energies (dict): A dictionary containing arrays of energies. The key 'total' should contain the total energy array.
    """
    if unit == "rydberg":
        factor = 1.0
        print("# potential energies (Ry)")
    elif unit == "hartree":
        factor = 0.5
        print("# potential energies (Ha)")
    elif unit == "eV":
        factor = Ry
        print("# potential energies (eV)")
    else:
        raise RuntimeError("Invalid unit: %s" % unit)
    print("#", end="")

    keys = list(energies.keys())
    key_list = ["total"] + list(keys[:-1])

    for i, key in enumerate(key_list):
        if i == 0:
            print("{:>16s}".format(key), end="")
        else:
            print("{:>17s}".format(key), end="")
    print("")

    energies_array = np.zeros((energies["total"].shape[0], len(key_list)), dtype=float)

    for i, key in enumerate(key_list):
        energies_array[:, i] = energies[key]

    for i in range(energies["total"].shape[0]):
        for j in range(len(key_list)):
            print("{:17.8e}".format(energies_array[i, j] * factor), end="")
        print("")


def print_forces(forces, unit="rydberg"):
    """
    Print the calculated forces.

    Parameters:
    - forces (numpy.ndarray): The array of forces.
    """

    if unit == "rydberg":
        factor = 1.0
        print("# atomic forces (Ry/Bohr)")
    elif unit == "hartree":
        factor = 0.5
        print("# atomic forces (Ha/Bohr)")
    elif unit == "eV":
        factor = Ry / Bohr
        print("# atomic forces (eV/Angstrom)")
    else:
        raise RuntimeError("Invalid unit: %s" % unit)

    forces_array = forces["total"] * factor
    for i in range(forces_array.shape[0]):
        print("# Snapshot %d" % (i + 1))
        for j in range(forces_array.shape[1]):
            print(
                "{:15.8e} {:15.8e} {:15.8e}".format(
                    forces_array[i, j, 0], forces_array[i, j, 1], forces_array[i, j, 2]
                )
            )


def read_dfset(fname, natom):
    """
    Read a DFSET file: displacements (Bohr), forces (Ry/Bohr), and, when every
    snapshot has an "E_pot (eV|Ry|Ha): value" comment, the DFT energies (Ry).
    """
    factors = {"(eV)": 1.0 / Ry, "(Ry)": 1.0, "(Ha)": 2.0}
    data, energies, pending = [], [], None
    with open(fname) as f:
        for line in f:
            if line.lstrip().startswith("#"):
                if "E_pot" not in line:
                    continue
                if pending is not None:
                    raise RuntimeError(
                        "Two E_pot headers for one snapshot in %s" % fname
                    )
                head, val = line.split("E_pot", 1)[1].split(":", 1)
                unit = [u for u in factors if u in head]
                if not unit:
                    raise RuntimeError("Unknown E_pot unit in %s" % line.strip())
                pending = float(val.split()[0]) * factors[unit[0]]
                if not np.isfinite(pending):
                    raise RuntimeError("Non-finite E_pot in %s" % line.strip())
            elif line.strip():
                data.append([float(t) for t in line.split()[:6]])
                if len(data) % natom == 0:  # end of a snapshot
                    energies.append(pending)
                    pending = None
    if len(data) % natom != 0:
        raise RuntimeError("The number of lines in %s is inconsistent." % fname)
    data = np.array(data).reshape(-1, natom, 6)
    # energies are compared only when every snapshot has one
    energies = None if None in energies else np.array(energies)
    return data[:, :, :3], data[:, :, 3:], energies


def compare_dfset(fname, energies, forces, e_dft, f_dft, e0, unit, plot):
    """
    Print the errors of the Taylor-expansion forces (and energies) against the
    DFT values of a DFSET file, and optionally save parity plots.
    """
    efac, elabel = {"rydberg": (1.0, "Ry"), "hartree": (0.5, "Ha"), "eV": (Ry, "eV")}[
        unit
    ]
    ffac, flabel = (
        efac / {"eV": Bohr}.get(unit, 1.0),
        elabel + ("/Angstrom" if unit == "eV" else "/Bohr"),
    )
    f_model = forces["total"]
    e_model = energies["total"]
    nsnap, natom = f_model.shape[:2]
    df = f_model - f_dft
    print("# %s: %d snapshots of %d atoms" % (fname, nsnap, natom))
    print(
        "# force error |F_model - F_DFT| / |F_DFT| = %.4f %%, RMS %.4e %s"
        % (
            100 * np.linalg.norm(df) / np.linalg.norm(f_dft),
            np.sqrt(np.mean(df**2)) * ffac,
            flabel,
        )
    )
    if e_dft is not None:
        # The DFT energy of the reference structure is not in the DFSET; unless
        # given by --e0, it is chosen so that the energy errors average to zero.
        e_ref = np.mean(e_dft - e_model) if e0 is None else e0 / Ry
        e_dft = e_dft - e_ref
        print(
            "# energy RMS error %.4e %s/atom (E0 = %.8f eV%s)"
            % (
                np.sqrt(np.mean((e_model - e_dft) ** 2)) / natom * efac,
                elabel,
                e_ref * Ry,
                "" if e0 is not None else ", fitted",
            )
        )
    print("#  snapshot  force error (%)", end="")
    if e_dft is not None:
        print("    E_DFT - E0 (%s)     E_model (%s)" % (elabel, elabel), end="")
    print("")
    for i in range(nsnap):
        print(
            "%11d %16.4f"
            % (i + 1, 100 * np.linalg.norm(df[i]) / np.linalg.norm(f_dft[i])),
            end="",
        )
        if e_dft is not None:
            print(" %17.8e %17.8e" % (e_dft[i] * efac, e_model[i] * efac), end="")
        print("")

    if plot is None:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    npanel = 1 if e_dft is None else 2
    fig, axes = plt.subplots(1, npanel, figsize=(4.5 * npanel, 4.2), squeeze=False)
    pairs = [(f_dft.ravel() * ffac, f_model.ravel() * ffac, "Force (%s)" % flabel)]
    if e_dft is not None:
        pairs.append((e_dft * efac, e_model * efac, "Energy (%s)" % elabel))
    for ax, (x, y, label) in zip(axes[0], pairs):
        lim = [min(x.min(), y.min()), max(x.max(), y.max())]
        ax.plot(lim, lim, color="0.6", lw=1)
        ax.plot(x, y, "o", ms=3, alpha=0.6)
        ax.set_xlabel("DFT " + label[0].lower() + label[1:])
        ax.set_ylabel("IFC model " + label[0].lower() + label[1:])
        ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(plot, dpi=150)
    print("# parity plot saved to %s" % plot)


def get_argoptions():
    """
    Parse and return command line arguments.

    Returns:
    - args: The parsed command line arguments.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--fcs", type=str, default=None, help="Force constant file")

    parser.add_argument(
        "--structure",
        metavar="SPOSCAR",
        default=None,
        help="The structure file that contains the information of original supercell "
        "(default: the supercell stored in the force constant file)",
    )

    parser.add_argument(
        "--format",
        metavar="vasp",
        default=None,
        help="The format of the structure specified by --structure.",
    )

    parser.add_argument(
        "--maxorder",
        type=int,
        default=4,
        help="The maximum order of the Taylor expansion. (2: Harmonic, 3: Cubic, etc.)",
    )

    parser.add_argument(
        "--disp",
        type=str,
        default=None,
        help="File containing the displacements in units of Bohr",
    )

    parser.add_argument(
        "--dfset",
        type=str,
        default=None,
        help="DFSET file (displacements in Bohr, forces in Ry/Bohr) whose DFT forces "
        "(and E_pot energies, if present) are compared with the IFC model",
    )

    parser.add_argument(
        "--e0",
        type=float,
        default=None,
        help="DFT energy (eV) of the undisplaced supercell for --dfset "
        "(default: chosen so that the energy errors average to zero)",
    )

    parser.add_argument(
        "--plot",
        metavar="FILE",
        default=None,
        help="With --dfset, save parity plots of the forces and energies to FILE "
        "(e.g. parity.png; requires matplotlib)",
    )

    parser.add_argument(
        "--unit",
        type=str,
        default="rydberg",
        choices=["rydberg", "hartree", "eV"],
        help="The unit of the output potential energies.",
    )

    parser.add_argument(
        "--print-forces",
        action="store_true",
        default=False,
        help="Print the forces in the output.",
    )

    return parser.parse_args()


def main():
    """
    Main function to execute the Taylor expansion potential calculation.
    """
    args = get_argoptions()
    parser = ForceConstantParser(args.fcs)
    if args.structure is None:
        lavec, xf, numbers = parser.supercell
        structure = Atoms(
            numbers=numbers, cell=lavec * Bohr, scaled_positions=xf, pbc=True
        )
    else:
        structure = get_original_structure(args.structure, args.format)

    fcs_dic = parser.read(maxorder=args.maxorder)
    taylor = TaylorExpansionPotential(structure, maxorder=args.maxorder)
    taylor.set_forceconstants(fcs_dic, parser.primitive_cell)

    natom = structure.get_global_number_of_atoms()
    if args.dfset is not None:
        displacements, f_dft, e_dft = read_dfset(args.dfset, natom)
        energies, forces = taylor.compute(displacements)
        compare_dfset(
            args.dfset, energies, forces, e_dft, f_dft, args.e0, args.unit, args.plot
        )
        return

    if args.disp is None:
        raise RuntimeError("Give --disp or --dfset.")

    displacements = np.loadtxt(args.disp)
    nlen = len(displacements) // (structure.get_global_number_of_atoms())
    if nlen * structure.get_global_number_of_atoms() != len(displacements):
        raise RuntimeError("The number of displacements is inconsistent.")

    displacements = displacements.reshape(-1, structure.get_global_number_of_atoms(), 3)
    energies, forces = taylor.compute(displacements)
    if args.print_forces:
        print_forces(forces, args.unit)
    else:
        print_energies(energies, args.unit)


if __name__ == "__main__":
    main()
