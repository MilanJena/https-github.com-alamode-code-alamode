This directory contains small python scripts and C++ programs which may be used as subsidiary tools.

## Python scripts

- displace.py : Generate input structure files of displaced configurations for VASP, Quantum-ESPRESSO, OpenMX, LAMMPS, and xTAPP.
- extract.py : Extract atomic displacements, forces, total energies, and effective charges from output files.
- plotband.py : Visualize phonon bands
- plotdos.py : Visualize phonon DOS and atom-projected phonon DOS
- analyzer.py : Compute phonon lifetimes, thermal conductivity, and cumulative thermal conductivity from a PREFIX.result file.
- scph_to_qefc.py : Create a new Quantum-ESPRESSO force constant file (q2r.x text or XML format) with the anharmonic correction of PREFIX.scph_dfc2 or PREFIX.scph.h5 (QHA likewise).
- dfc2.py : Apply SCPH/QHA harmonic corrections (PREFIX.scph_dfc2 / PREFIX.qha_dfc2) to an ALAMODE HDF5 force-constant file, or extract the renormalized FC2 at one temperature from PREFIX.scph.h5 / PREFIX.qha.h5.
- fcsinfo.py : Print the content of an ALAMODE force-constant HDF5 file or of PREFIX.scph.h5 / PREFIX.qha.h5 (structure, FC orders, distances, ASR residual, provenance, SCPH/QHA settings).
- fcconvert.py : Convert harmonic and cubic force constants between ALAMODE (.h5/.xml), phonopy/phono3py, and ShengBTE formats, and quartic force constants between ALAMODE and FourPhonon. FC3 periodic images are exact only for alamode -> ShengBTE; alamode -> phonopy/phono3py is exact only at commensurate q when FC3 reaches the supercell boundary (see the docstring).
- calcpes.py, taylor.py : Evaluate the Taylor-expansion potential (energies and forces) defined by ALAMODE force constants, and compare it with DFT forces and energies of a DFSET file (--dfset, parity plot with --plot).
- strainifc.py : Generate the strain-IFC coupling inputs (``strain_harmonic.in``, ``strain_force.in``) for the SCPH/QHA structural optimization with cell relaxation in anphon. Port of the strainIFCcoupling scripts by Ryota Masuki (https://github.com/r-masuki/strainIFCcoupling).
- elastic.py : DFT finite-strain workflow for the first-, second-, and third-order elastic constants (``elastic_constants.in``, ``C1_array.in``) used by the same optimization.

The ``strainkit/`` package contains the shared code of ``strainifc.py`` and ``elastic.py``; ``fcsio/``, ``interface/``, and ``analyzer/`` are helper packages of the other scripts.

To use the scripts, Python environment (+ Numpy) is necessary.
Matplotlib is also required for plotband.py and plotdos.py; ase and spglib for strainifc.py, elastic.py, calcpes.py, and the fcsio package (which also needs h5py and lxml, as does fcconvert.py);
the ``alm`` Python package built from the ``python/`` directory of ALAMODE for the harmonic fits of strainifc.py.

To see available options of each script, please run the script with ``--help`` option.

Tests of the strainkit package: ``cd tools && python -m pytest tests``.

## C++ scripts (binaries)

- analyze\_phonons.\* : Computes (and prints out) phonon lifetimes, mean-free-path, and (cumulative) thermal-conductivity using the .result file as an input
- qe2alm.\* : Converts a Quantum-ESPRESSO force constant file to the ALAMODE XML format.
- dfc2.\* : Create effective harmonic force constant (ALAMODE XML) from input harmonic force constant (XML format) and PREFIX.scph_dfc2.
- virtual.\* : Performs linear interpolation of force constants

To use these code, please edit the Makefile and do make (or use Cmake).
