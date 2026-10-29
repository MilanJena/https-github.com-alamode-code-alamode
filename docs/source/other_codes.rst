.. _label_other_codes:

Using other codes (ASE calculators, WIEN2k, CRYSTAL, ...)
=========================================================

.. |Angstrom|   unicode:: U+00C5

``displace.py`` and ``extract.py`` support VASP, Quantum ESPRESSO, OpenMX, xTAPP, and LAMMPS.
Any other code that calculates atomic forces can be used as well.
*alm* only needs a ``DFSET`` file: the displacements and forces of each structure, in Cartesian coordinates
(see :ref:`the DFSET format <label_format_DFSET>`).
This page shows how to make the structures and the ``DFSET`` file
with `ASE <https://wiki.fysik.dtu.dk/ase/>`_, and what to look out for with codes such as WIEN2k and CRYSTAL.


The general recipe
------------------

1. **Make the displaced structures in the VASP format**, whatever code you use.
   The POSCAR files are only an exchange format: they can be read by ASE and converted to the input of other codes.
   With :ref:`SUPERCELL <alm_supercell>`, give *alm* the unit cell and run::

     $ alm alm0.in > alm0.log          # MODE = suggest
     $ displace.py --VASP POSCAR --mag 0.01 --prefix disp --pattern_file PREFIX.pattern_HARMONIC

   This writes the displaced supercells ``disp*.POSCAR`` and the undisplaced supercell ``SPOSCAR0``.
   If you build the supercell yourself, give it to *alm* and to ``displace.py`` instead; it is then the reference structure in place of ``SPOSCAR0``.
   For anharmonic force constants, ``displace.py --VASP SPOSCAR0 --random_normalcoord ...`` (or ``--random``) writes the structures in the same way.

2. **Calculate the forces** of each structure with your code (sections below).
   These are single-point calculations: do not relax the atomic positions or the cell.

3. **Write the DFSET file.**
   Take the displacements from the POSCAR files (displaced structure minus ``SPOSCAR0``), not from the output of your code,
   and take only the forces from your code.
   The atoms must be in the order of ``SPOSCAR0``, and the forces must be in the Cartesian axes of ``SPOSCAR0``.

4. Run *alm* with ``MODE = optimize`` and ``DFSET`` as usual.

The rest of the workflow (*anphon*, the analysis scripts) does not depend on the force calculator.
For polar materials, the dielectric constant and the Born effective charges from your code are given to *anphon* in the
:ref:`BORNINFO format <label_format_BORNINFO>`.


ASE calculators
---------------

ASE provides calculators for many codes, for example GPAW, ABINIT, CP2K, CASTEP, FHI-aims, SIESTA, and CRYSTAL,
as well as machine-learning potentials.
The script below calculates the forces of all ``disp*.POSCAR`` files and writes the ``DFSET`` file in bohr and Ry/bohr,
the default units of *alm*.
Replace ``EMT()`` with the calculator of your choice::

  import glob

  import numpy as np
  from ase.calculators.emt import EMT  # replace with your calculator
  from ase.io import read
  from ase.units import Bohr, Rydberg


  def get_forces(atoms):
      atoms.calc = EMT()
      return atoms.get_forces()  # eV/Angstrom


  ref = read("SPOSCAR0", format="vasp")
  f0 = get_forces(ref.copy())  # residual forces of the undisplaced supercell

  with open("DFSET", "w") as out:
      for name in sorted(glob.glob("disp*.POSCAR")):
          atoms = read(name, format="vasp")
          forces = get_forces(atoms) - f0
          d = atoms.get_scaled_positions(wrap=False) - ref.get_scaled_positions(wrap=False)
          d -= np.round(d)  # nearest periodic image
          u = d @ ref.cell.array  # Angstrom
          out.write(f"# {name}\n")
          for ui, fi in zip(u / Bohr, forces / (Rydberg / Bohr)):
              out.write("%20.14f %20.14f %20.14f %20.12e %20.12e %20.12e\n" % (*ui, *fi))

Subtracting the forces of the undisplaced supercell (``f0``) removes small residual forces, as ``extract.py --offset`` does.
It can be dropped if the structure is well relaxed.
``ase.units`` uses CODATA 2014 by default, whereas ALAMODE uses CODATA 2018; the difference (about 10\ :sup:`-8`) is negligible.
Alternatively, write the data in |Angstrom| and eV/|Angstrom| and set
``LENGTH_UNIT = angstrom`` and ``FORCE_UNIT = eV/angstrom`` in *alm*
(note that :ref:`LENGTH_UNIT <alm_length_unit>` also applies to ``&cell`` and ``&cutoff``).

As a check, the script above with the EMT potential of ASE for fcc Cu (``SUPERCELL = 3`` of the primitive cell, 27 atoms)
gives a fitting error of 0.9 %, and the phonon frequencies at X and L agree with those of ``ase.phonons`` to within 0.5 %.

The force constants can also be fitted without files, through the :doc:`alm Python package <alm_python>`.

.. note::

   **CRYSTAL** is easiest to use through ``ase.calculators.crystal``.
   It writes the structure to ``fort.34`` (``EXTERNAL`` geometry) with the full lattice vectors, in the given atom order,
   and with the identity as the only symmetry operation, and calculates the forces with ``GRADCAL``.
   The forces are converted to eV/|Angstrom| by ASE.
   See the ASE documentation for the basis-set file and the command to run CRYSTAL.


Codes without an ASE calculator (e.g. WIEN2k)
---------------------------------------------

For other codes, convert the POSCAR files to the input format of the code
(``ase.io.write`` supports many formats, or use the tools of the code), run the calculations,
read the forces from the output, and write the ``DFSET`` file with the same loop as above.
Check the following points, because a mistake here does not stop *alm*; it only gives wrong force constants.

* **Atom order.** Some codes reorder the atoms, for example by grouping atoms that are equivalent by symmetry.
  Turn this off, or map the atoms back by their positions, comparing the structure read from the output of the code
  with the displaced POSCAR it was made from.
  Fractional coordinates do not depend on the Cartesian axes, so the mapping can be done with them,
  provided that the code keeps the lattice vectors (in the same order) and the origin::

    def order_map(ref, atoms, tol=1.0e-3):
        """atoms[order_map(ref, atoms)[i]] is atom i of ref."""
        d = atoms.get_scaled_positions()[None, :, :] - ref.get_scaled_positions()[:, None, :]
        d -= np.round(d)
        dist = np.linalg.norm(d @ ref.cell.array, axis=2)
        order = np.argmin(dist, axis=1)
        assert len(set(order)) == len(ref), "the mapping is not one-to-one"
        assert dist[np.arange(len(ref)), order].max() < tol, "the positions do not match"
        assert (atoms.numbers[order] == ref.numbers).all(), "the species do not match"
        return order

    forces = forces[order_map(atoms, atoms_from_output)]  # atoms: the displaced POSCAR

* **Cartesian axes.** Codes that build the lattice from the cell parameters :math:`(a, b, c, \alpha, \beta, \gamma)`
  use their own Cartesian axes, which can differ from those of ``SPOSCAR0`` by a rotation.
  This matters for non-orthogonal supercells, such as a supercell of an fcc primitive cell.
  Rotate the forces back with the lattice vectors of the code ``A_code`` (rows are the lattice vectors, in the axes of the code,
  in the same order as in ``SPOSCAR0``)::

    A_ref = ref.cell.array
    R = A_ref.T @ np.linalg.inv(A_code.T)
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-4)  # otherwise the cell has changed
    forces = forces @ R.T

  Alternatively, use a supercell whose lattice vectors are along *x*, *y*, and *z* (e.g. a diagonal ``SUPERCELL`` of a cubic or tetragonal cell).

* **Units.** Convert the forces to Ry/bohr (or declare the unit with :ref:`FORCE_UNIT <alm_force_unit>`).

* **Symmetry.** A displaced supercell has lower symmetry than the perfect crystal.
  If the code symmetrizes the forces or the charge density, it must use only the symmetry of the displaced structure.

.. note::

   **WIEN2k.** ``ase.io.write("case.struct", atoms, format="struct")`` writes every atom as an inequivalent atom (``MULT = 1``)
   with lattice type P. Set the muffin-tin radii (the ``rmt`` argument, a list with one value in bohr per atom, or ``setrmt_lapw``) so that the spheres do not touch after the displacement.
   Keep the atoms inequivalent and in this order when initializing the calculation (do not accept a structure from ``sgroup`` that merges equivalent atoms);
   check that ``case.scf`` lists the forces of all atoms.
   Turn on the force calculation (e.g. ``run_lapw -fc ...``) and read the total forces of the last iteration
   from the ``:FGL`` lines of ``case.scf``. Check their unit (mRy/bohr) and Cartesian axes in the WIEN2k user guide for your version;
   WIEN2k builds the lattice from the cell parameters, so rotate the forces as above for non-orthogonal cells.
