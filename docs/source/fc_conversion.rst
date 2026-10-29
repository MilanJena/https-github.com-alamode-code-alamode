.. _label_fc_conversion:

Converting force constants (phonopy, phono3py, ShengBTE)
========================================================

.. |Angstrom|   unicode:: U+00C5

``tools/fcconvert.py`` converts harmonic (FC2) and cubic (FC3) force constants between ALAMODE,
phonopy/phono3py, and ShengBTE, in any direction, and quartic (FC4) force constants between ALAMODE and FourPhonon,
the four-phonon extension of ShengBTE.
Typical uses are calculating the thermal conductivity of ALAMODE force constants with phono3py or ShengBTE,
or using force constants of phonopy or phono3py in *anphon* (``FCSFILE``, ``FC2FILE``, ``FC3FILE``).

The script needs Python 3 with numpy, h5py, ASE, spglib, and lxml (the last two through ``tools/fcsio``),
and PyYAML to read the phonopy and phono3py yaml files.
Run it as ``tools/fcconvert.py`` (or put ``tools/`` on your ``PATH``) so that the ``fcsio`` package next to it is found.
``fcconvert.py --help`` prints the full description.


Supported files
---------------

``--from`` and ``--to`` take ``alamode``, ``phonopy``, ``phono3py``, or ``shengbte``; ``phono3py`` is an alias of ``phonopy``.
The force constants are given with ``--fc2``, ``--fc3``, and/or ``--fc4`` (``alamode`` and ``shengbte`` only), and the output directory with ``-o`` (default: the current directory).

.. list-table::
   :header-rows: 1
   :widths: 14, 43, 43

   * - Format
     - Read
     - Written
   * - ``alamode``
     - ``PREFIX.h5`` or legacy ``PREFIX.xml`` of *alm*, *anphon*, or this script.
       No ``--structure`` needed.
     - ``PREFIX.h5`` (``--prefix``, default ``fcconvert``), readable by *anphon* as ``FCSFILE``/``FC2FILE``/``FC3FILE``.
       If the orders have different supercells, one file per order (``PREFIX_fc2.h5``, ``PREFIX_fc3.h5``, ``PREFIX_fc4.h5``).
   * - ``phonopy`` / ``phono3py``
     - FC2: ``FORCE_CONSTANTS``, ``force_constants.hdf5``, or ``fc2.hdf5``; FC3: ``fc3.hdf5``; full or compact.
       ``--structure``: ``phonopy.yaml``, ``phonopy_disp.yaml``, ``phono3py.yaml``, ``phono3py_disp.yaml``,
       or a POSCAR-type supercell (``SPOSCAR``) together with ``--supercell``.
     - With ``--fc2``: ``FORCE_CONSTANTS`` (``force_constants.hdf5`` with ``--hdf5``) and ``phonopy.yaml``.
       With ``--fc3``: in addition ``fc2.hdf5`` (if FC2 is given), ``fc3.hdf5``, and ``phono3py.yaml``.
       Full arrays by default, compact ones with ``--compact``.
   * - ``shengbte``
     - ``FORCE_CONSTANTS_2ND``, ``FORCE_CONSTANTS_3RD``, ``FORCE_CONSTANTS_4TH`` (FourPhonon); ``--structure``: ``CONTROL``.
     - ``FORCE_CONSTANTS_2ND`` (with ``--fc2``), ``FORCE_CONSTANTS_3RD`` (with ``--fc3``), ``FORCE_CONSTANTS_4TH`` (with ``--fc4``),
       and a minimal ``CONTROL``.

``--supercell`` takes 3 (diagonal) or 9 integers in the convention of phonopy's ``DIM``.
It is required when reading a POSCAR-type supercell, and when writing ``FORCE_CONSTANTS_2ND`` from a supercell that is not
diagonal in the primitive cell (see the first example). With ``--fc3`` alone, ``scell`` in the written ``CONTROL`` is
taken from ``--supercell`` and is ``1 1 1`` without it, so set it to the supercell of your harmonic force constants.


Examples
--------

The examples use the harmonic and cubic force constants of Si in the 2\ |times|\ 2\ |times|\ 2 conventional supercell
(64 atoms) of the :ref:`silicon tutorial <label_tutorial_01>`, ``si222.h5`` and ``si222_cubic.h5``.

.. |times| unicode:: U+00D7

**ALAMODE to ShengBTE** ::

  $ fcconvert.py --from alamode --to shengbte --fc2 si222.h5 --fc3 si222_cubic.h5 --supercell 4 4 4 -o shengbte

ShengBTE needs a supercell of :math:`n_1\times n_2\times n_3` primitive cells, which the conventional supercell is not.
``--supercell 4 4 4`` refolds the harmonic force constants onto a 4\ |times|\ 4\ |times|\ 4 primitive supercell (128 atoms).
Choose a supercell that is at least as large as the range of the force constants.
The written ``CONTROL`` holds only the crystal structure, ``scell``, and placeholder values
(``ngrid = 10 10 10``, ``T = 300``, ``nonanalytic = .FALSE.``): edit it before running ShengBTE,
and add the dielectric tensor and the Born effective charges for a polar material.
For FourPhonon, add ``--fc4`` with the file of the quartic force constants (the same file as ``--fc3`` if it holds both)
and set ``four_phonon = .TRUE.`` in ``CONTROL``.

**ALAMODE to phono3py** ::

  $ fcconvert.py --from alamode --to phono3py --fc2 si222.h5 --fc3 si222_cubic.h5 -o phono3py

This writes ``fc2.hdf5``, ``fc3.hdf5``, and ``phono3py.yaml`` (and ``FORCE_CONSTANTS`` and ``phonopy.yaml`` for phonopy).
The yaml files take the primitive cell as the unit cell (``primitive_matrix`` is the identity),
and the supercell matrix is that of the ALAMODE supercell.
With ``--to phonopy --fc2 si222.h5`` only the harmonic files are written.

**phonopy or phono3py to ALAMODE** ::

  $ fcconvert.py --from phonopy --to alamode --fc2 FORCE_CONSTANTS --structure phonopy_disp.yaml --prefix si_phonopy
  $ fcconvert.py --from phono3py --to alamode --fc2 fc2.hdf5 --fc3 fc3.hdf5 --structure phono3py_disp.yaml --prefix si_phono3py

The first command writes ``si_phonopy.h5``, which *anphon* reads with ``FCSFILE = si_phonopy.h5`` (or ``FC2FILE``).
The primitive cell is taken from the yaml file (``primitive_cell`` or ``primitive_matrix``).
For phono3py, the harmonic force constants belong to ``phonon_supercell`` when the yaml file defines one.
If it differs from the FC3 supercell, FC2 and FC3 go into separate files ``si_phono3py_fc2.h5`` and ``si_phono3py_fc3.h5``,
to be given to *anphon* as ``FC2FILE`` and ``FCSFILE``.
Only yaml files in |Angstrom| and eV/|Angstrom|\ :sup:`2` (VASP-like units) are accepted.
With a supercell POSCAR instead of a yaml file, give the ``DIM`` with which phonopy made it::

  $ fcconvert.py --from phonopy --to alamode --fc2 FORCE_CONSTANTS --structure SPOSCAR --supercell 2 2 2 --prefix si_phonopy

The unit cell of the phonopy calculation is then used as the primitive cell. If that is a conventional cell,
prefer the yaml file, which also gives the primitive cell.

**ShengBTE to ALAMODE** ::

  $ fcconvert.py --from shengbte --to alamode --fc2 FORCE_CONSTANTS_2ND --fc3 FORCE_CONSTANTS_3RD --structure CONTROL

``FORCE_CONSTANTS_3RD`` carries no supercell; it is assigned the ``scell`` of ``CONTROL`` (or ``--supercell``).
Add ``--fc4 FORCE_CONSTANTS_4TH`` to include the quartic force constants of FourPhonon.


Conventions
-----------

* **Units.** ALAMODE ``.xml`` files are in Ry/bohr\ :sup:`n`; ``.h5`` files are read according to their unit attributes.
  Written ``.h5`` files are in eV/|Angstrom|\ :sup:`n` with unit attributes, which *anphon* converts.
  phonopy, phono3py, and ShengBTE use eV/|Angstrom|\ :sup:`n`; the lattice vectors in ``FORCE_CONSTANTS_3RD`` and ``FORCE_CONSTANTS_4TH`` are in |Angstrom|,
  and ``CONTROL`` holds the same |Angstrom| values with ``lfactor = 0.1`` (ShengBTE multiplies them by ``lfactor`` to get nm),
  so keep ``lfactor`` if you edit the vectors.
* **Atom order.** Written supercells follow phonopy's default atom order: atom-major, with the lattice points running fastest along :math:`x`.
  This is also the order ShengBTE expects for ``FORCE_CONSTANTS_2ND``, but it generally differs from the order of the ALAMODE supercell.
  Use the written files together with the ``phonopy.yaml``, ``phono3py.yaml``, or ``CONTROL`` written next to them,
  not with the ``SPOSCAR`` or the *alm* input of the ALAMODE calculation.
  When you read phonopy files, the ``--structure`` file must be the one that belongs to the force constants.
* **Primitive cell.** Force constants given together must use the same primitive cell with the same atom order;
  otherwise the script stops.


Periodic images
---------------

A supercell tensor (phonopy, phono3py, ``FORCE_CONSTANTS_2ND``) assigns one value to each pair or triplet of supercell atoms;
it does not say which periodic image of an atom interacts.
ALAMODE, ``FORCE_CONSTANTS_3RD``, and ``FORCE_CONSTANTS_4TH`` list the periodic images explicitly.
When ``fcconvert.py`` reads a supercell tensor, it assigns each pair all the minimum-distance images (within 10\ :sup:`-4` |Angstrom|)
with the value split equally; for FC3 the images of the second and third atoms are chosen independently.
This reproduces the supercell tensor, and hence the dynamical matrix and the anharmonic couplings at every :math:`q` point
commensurate with the supercell, exactly. Beyond that:

* ALAMODE to ShengBTE keeps the images of FC3 and FC4 exactly, and ALAMODE to ALAMODE reproduces the original images.
* Harmonic force constants are unaffected in practice, because ALAMODE also splits FC2 over the minimum-distance images.
  For the Si example, an ALAMODE |rarr| phonopy |rarr| ALAMODE round trip of FC2 gives the same *anphon* frequencies,
  also at non-commensurate :math:`q`.
* ALAMODE assigns the images of a cubic cluster jointly, which a supercell tensor cannot represent.
  ALAMODE to phonopy/phono3py is therefore exact only at :math:`q` points commensurate with the supercell when FC3 reaches the supercell boundary
  (long or no cutoff). The thermal conductivity on non-commensurate meshes then changes by about 1% at 300 K
  and 3-5% at 100 K (Si 2\ |times|\ 2\ |times|\ 2 conventional cell, no FC3 cutoff), and phono3py's default
  ``r0`` averaging adds another 1-1.5% (its ``--v2`` option reproduces the convention described above).
  FC3 with a cutoff shorter than half the supercell is unaffected.
* Force constants of phonopy, phono3py, or ShengBTE read into ALAMODE follow the image convention of the source code, not that of ALAMODE.

.. |rarr| unicode:: U+2192


.. _label_fc_conversion_alm:

Writing ShengBTE files directly from *alm*
------------------------------------------

*alm* can write the anharmonic force constants of a fit in the ShengBTE format itself:
:ref:`FC3_SHENGBTE = 1 <alm_fc3_shengbte>` writes ``PREFIX.FORCE_CONSTANTS_3RD`` (``NORDER >= 2``), and
:ref:`FC4_SHENGBTE = 1 <alm_fc4_shengbte>` writes ``PREFIX.FORCE_CONSTANTS_4TH`` (``NORDER >= 3``) for FourPhonon,
the four-phonon extension of ShengBTE.
Remove the ``PREFIX.`` part of the names (or link the files) before running ShengBTE or FourPhonon.
These files have the same image convention as ``fcconvert.py --to shengbte``, but *alm* writes neither
``FORCE_CONSTANTS_2ND`` nor ``CONTROL``; use ``fcconvert.py`` for those, or to convert force constants you already have
(``--fc4`` for ``FORCE_CONSTANTS_4TH``).
