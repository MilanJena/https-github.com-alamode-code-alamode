.. _label_faq:

Frequently Asked Questions (FAQ)
================================

This page collects questions that come up repeatedly on the
`GitHub Discussions <https://github.com/alamode-code/alamode/discussions>`_ and
`issues <https://github.com/alamode-code/alamode/issues>`_ pages, with answers for
ALAMODE 2.0. Click a question to show the answer. If your question is not answered
here, search the Discussions first, then open a new thread with the ALAMODE version,
the compiler and MPI library, the input files, and the error message.

.. grid:: 1 2 3 3
   :gutter: 2

   .. grid-item-card:: :octicon:`tools;1.2em` General & installation
      :link: faq_general
      :link-type: ref

      Building, crashes at startup, units and conventions.

   .. grid-item-card:: :octicon:`graph;1.2em` Force constants (alm)
      :link: faq_alm
      :link-type: ref

      Fitting error, cutoffs, supercells, cross-validation, symmetry.

   .. grid-item-card:: :octicon:`pulse;1.2em` Harmonic phonons
      :link: faq_phonons
      :link-type: ref

      Imaginary modes, LO-TO splitting, 2D materials, band plots.

   .. grid-item-card:: :octicon:`flame;1.2em` Thermal conductivity
      :link: faq_kappa
      :link-type: ref

      Convergence, BTE solvers, four-phonon, lifetimes.

   .. grid-item-card:: :octicon:`sun;1.2em` Finite temperature
      :link: faq_finite_t
      :link-type: ref

      SCPH convergence, imaginary modes, renormalized force constants.

   .. grid-item-card:: :octicon:`arrow-switch;1.2em` Interfaces & performance
      :link: faq_interfaces
      :link-type: ref

      Other codes, file conversion, memory, MPI and OpenMP.


.. _faq_general:

General and installation
------------------------

.. dropdown:: The build fails because spglib (or another library) is not found.
   :icon: question

   The easiest route is the conda or pixi environment described in :doc:`install`, which provides all
   dependencies (compilers, MPI, Boost, Eigen, spglib, HDF5, ...) in one place; CMake then finds them in
   ``$CONDA_PREFIX`` automatically. Outside such an environment, or if CMake still cannot find spglib
   (e.g. ``cannot find -lsymspg``), give its location explicitly::

      $ cmake -DSPGLIB_ROOT=/path/to/spglib -DHDF5_ROOT=/path/to/hdf5 ..

   ``SPGLIB_ROOT`` is the directory that contains the ``lib`` and ``include`` directories of spglib.
   If the binaries cannot find the shared libraries at run time, add them to the library path::

      $ export LD_LIBRARY_PATH=$SPGLIB_ROOT/lib:$HDF5_ROOT/lib:$LD_LIBRARY_PATH

   See also `issue #130 <https://github.com/alamode-code/alamode/issues/130>`_.

.. dropdown:: The program stops with "Illegal instruction" or crashes in the MPI initialization.
   :icon: question

   * **Illegal instruction.** A build with ``-xHOST`` (Intel) or ``-march=native`` (GCC) uses the instruction
     set of the machine where it was compiled. Running it on an older CPU, e.g. on other nodes of a cluster,
     stops with this error. Build on the machine type where the code runs, or use a more generic option.
   * **Crash in MPI_Init / MPI_Comm_rank.** This usually means that the ``mpirun`` command and the MPI library
     linked to the binary come from different MPI installations (e.g. a system MPI and the conda OpenMPI).
     Use the ``mpirun`` of the MPI library used for the build, and check the linked libraries with ``ldd``.

.. dropdown:: anphon spends several seconds at startup before any output when many MPI processes run on one node. Why?
   :icon: question

   This comes from the OpenMP runtime of the Intel oneAPI 2025 compilers (``libiomp5``), not from **anphon**.
   When each MPI process enters its first OpenMP parallel region, the runtime reads the hardware topology of
   the node through the hwloc library. When dozens of processes do this at the same time, the start of the run
   is delayed by a few seconds (about 3--6 seconds with 64 processes on one node), independently of the input.
   The OpenMP runtimes of GCC (``libgomp``), LLVM (``libomp``), and Intel oneAPI 2023 do not show this delay.

   The delay matters only for short runs. If you want to avoid it with the oneAPI 2025 runtime, let hwloc read
   the topology from a file generated on the compute node at the start of the job script::

      $ lstopo-no-graphics --of xml > topo_${JOB_ID}.xml
      $ export HWLOC_XMLFILE=${PWD}/topo_${JOB_ID}.xml HWLOC_THISSYSTEM=1
      $ mpirun -np 64 anphon input.in > output.log

   ``lstopo-no-graphics`` is part of hwloc. Generate the file on the node where **anphon** runs, because a
   topology file from a different type of node gives wrong information to the runtime.

.. dropdown:: Which units and coordinate frames are used?
   :icon: question

   * *alm* works in Rydberg atomic units: lattice vectors and displacements in bohr, forces in Ry/bohr, and
     :math:`n`-th order force constants in Ry/bohr\ :sup:`n`. Data in other units can be given with
     :ref:`LENGTH_UNIT <alm_length_unit>` and :ref:`FORCE_UNIT <alm_force_unit>`, e.g.
     ``LENGTH_UNIT = angstrom; FORCE_UNIT = eV/angstrom`` for a ``DFSET`` in VASP units.
     ``extract.py`` writes Rydberg atomic units by default.
   * Vectors and tensors in the outputs (force constants in ``PREFIX``.h5, eigenvectors, group velocities,
     mean square displacements, thermal conductivity) are given in the Cartesian frame, not along the lattice
     vectors.
   * Phonon frequencies are given in cm\ :sup:`-1`. ``plotband.py --unit`` converts the plotted bands to other
     units (see ``plotband.py -h``).
   * Unit conversions use the CODATA 2018 values of the physical constants, in the C++ codes
     (``include/constants.h``) and in the Python tools (``tools/alamode_constants.py``) alike.


.. _faq_alm:

Force constants with alm
------------------------

.. dropdown:: The fitting error is very large (> 90%). Is it problematic? If so, how can I reduce the error?
   :icon: question

   A large fitting error can affect the accuracy of the force constants. The most likely cause is non-zero
   residual forces in the original supercell structure (before the atoms are displaced). Even when the
   structure of the primitive cell is optimized with a relatively strict convergence criterion, the atomic
   forces in the supercell generated from it may deviate from zero. To reduce the error associated with the
   residual forces, calculate the forces of the undisplaced supercell as well and subtract them with the
   ``--offset`` option of ``extract.py``. For VASP, for example::

      $ extract.py --VASP SPOSCAR0 --offset vasprun0.xml vasprun_harm*.xml > DFSET_harmonic

   Here, ``vasprun0.xml`` is the output of the VASP calculation for the undisplaced supercell ``SPOSCAR0``.

   If the fitting error is still large after subtracting the offset components, please consider the
   following points:

   1. Give the lattice vectors and fractional coordinates with ~15 significant digits; for example, 1/3
      should be 0.333333333333333 rather than 0.33333. Alternatively, let *alm* symmetrize the structure
      with :ref:`SYMMETRIZE_STRUCTURE = 1 <alm_symmetrize_structure>`.
   2. Check that the order of the atoms in the ``DFSET`` file is the same as in the supercell of *alm*.
      This matters when the supercell or the displacements were made by another program.
   3. Check that the DFT calculations are properly converged.
   4. Use a smaller displacement magnitude for the harmonic force constants.

.. dropdown:: How small should the fitting error be?
   :icon: question

   It depends on the Taylor expansion potential and displacement magnitude you choose.

   In a standard harmonic calculation where ``--mag 0.01`` is used in ``displace.py`` and all harmonic
   interactions are considered (cutoff = None), the fitting error is usually less than 5% (~1--2% in most
   cases).

   In the calculation of third-order force constants with ``--mag 0.04``, the fitting error should be small
   as well. Indeed, in many cases, the fitting errors are much smaller (< 1%) than in the harmonic case.

   When the force constants are fitted to snapshots of molecular dynamics or random displacements at finite
   temperature, the fitting error is larger (often > 10%). In this case, compare the training and validation
   errors:

   * Both errors are large: the model is not flexible enough. Increase the cutoff radii or ``NORDER``.
   * The training error is small but the validation error is much larger: the model is overfitted.
     Add more training data.

   See `discussion #258 <https://github.com/alamode-code/alamode/discussions/258>`_ and
   `discussion #328 <https://github.com/alamode-code/alamode/discussions/328>`_.

.. dropdown:: What value should I use for the cutoff radius?
   :icon: question

   For the harmonic term, we recommend ``None``, which includes all harmonic interactions inside the
   supercell. This choice hardly increases the computational cost because the number of displacement
   patterns does not change. Also, the harmonic dynamical matrix becomes exact at the commensurate q points
   only when ``None`` is selected. (Giving a very large cutoff radius has the same effect as giving ``None``.)

   For the anharmonic terms, increase the cutoff radii gradually and check the convergence of physical
   quantities, such as thermal conductivity and free energy, with respect to the cutoff values. In most
   cases, a cutoff radius of 10 Bohr is a good starting point, but a larger value may be necessary for polar
   materials. The cutoffs are not chosen automatically.

   If *alm* prints ``WARNING : Rank deficient``, the training data are not enough to determine all
   parameters of the model. Reduce the cutoff radii (or :ref:`NBODY <alm_nbody>`), or add more data,
   e.g. the structures suggested by ``MODE = suggest`` for the present cutoffs, the same patterns with a
   different ``--mag``, or more random displacements.

.. dropdown:: How do I choose NORDER and NBODY?
   :icon: question

   :ref:`NORDER <alm_norder>` = :math:`m` includes the terms up to the :math:`(m+1)`-th order: ``NORDER = 2``
   (up to cubic) is enough for three-phonon thermal conductivity, and ``NORDER = 3`` (up to quartic) is needed
   for SCPH and four-phonon scattering. A higher ``NORDER`` can be used as a more
   flexible model when the fitting error stays large for strongly anharmonic materials or high-temperature
   data.

   :ref:`NBODY <alm_nbody>` removes many-body clusters, whose number grows quickly. ``NORDER = 3; NBODY = 2 3 3``
   (no four-body quartic terms) is a common choice. Harmonic force constants always involve two atoms, so the
   first entry of ``NBODY`` is always 2.

.. dropdown:: How large should the displacements be?
   :icon: question

   * Harmonic force constants: about 0.01 Å (``displace.py --mag 0.01``). Very small displacements
     (e.g. 10\ :sup:`-5` Å) make the forces dominated by numerical noise of the DFT calculation.
   * Anharmonic force constants from systematic displacements: 0.03--0.05 Å usually works well, but check how
     the value affects your results.
   * Random displacements: use ``displace.py --VASP SPOSCAR0 --prim POSCAR --random_normalcoord --evec PREFIX.evec.h5 --temp T``
     (the primitive cell given with ``--prim`` is required), which samples the displacements from the harmonic phonons at temperature ``T`` (see :doc:`quickstart` and
     :ref:`label_tutorial_silicon_ifc`).

   See `discussion #264 <https://github.com/alamode-code/alamode/discussions/264>`_.

.. dropdown:: Can I use different supercells for the harmonic and anharmonic force constants?
   :icon: question

   Yes. Anharmonic interactions are usually shorter ranged than harmonic ones, so the anharmonic force
   constants are often calculated in a smaller supercell.

   * In *anphon*, give the harmonic force constants of the larger supercell with
     :ref:`FC2FILE <anphon_fc2file>` and the anharmonic ones with ``FCSFILE`` (or ``FC3FILE``/``FC4FILE``).
   * In *alm*, the file given to :ref:`FC2FIX <alm_fc2fix>` must contain harmonic force constants of the
     **same** supercell as the present fit. To fix the harmonic part while fitting the anharmonic terms in a
     small supercell, calculate the harmonic force constants of that small supercell first.

   Check the convergence with respect to the supercell size, in particular along short lattice directions:
   a supercell that is too short in one direction can underestimate the phonon-phonon scattering.

.. dropdown:: How do I use the elastic net or adaptive LASSO, and how do I choose L1_ALPHA?
   :icon: question

   1. Run cross-validation first: ``LMODEL = enet`` (or ``adaptive-lasso``) with ``CV = 4`` or so
      (automatic ``CV``-fold cross-validation), or ``CV = -1`` with a separate ``DFSET_CV`` (manual).
   2. Take the ``L1_ALPHA`` that minimizes the validation error in ``PREFIX``.cvscore (or ``PREFIX``.cvset).
      If the minimum is at the smallest alpha, rerun with a smaller :ref:`CV_MINALPHA <alm_cv_minalpha>`.
   3. Run the fit with ``CV = 0`` and that ``L1_ALPHA``. The force-constant files are written only by this
      run.

   The ``RESIDUAL`` printed in the log is the training error. A smaller alpha always lowers it, so it cannot
   be used to choose alpha. The optimal alpha changes with the model (cutoffs, ``NORDER``, ``NBODY``,
   ``FCSYM_BASIS``), with the amount of training data, and with ``STANDARDIZE``, so repeat the
   cross-validation when any of them changes. Skipping the cross-validation with a tiny ``L1_ALPHA`` is not
   recommended because it overfits when the training data are few. ``DEBIAS_OLS = 1`` refits the non-zero
   coefficients by least squares.

   See :ref:`LMODEL <alm_lmodel>`, :ref:`CV <alm_cv>`, :ref:`label_tutorial_bto_ifc`, and
   `discussion #328 <https://github.com/alamode-code/alamode/discussions/328>`_.

.. dropdown:: alm takes very long (or crashes) while generating the constraints.
   :icon: question

   The most common cause is a structure that is symmetric only approximately, e.g. coordinates given with
   too few digits (hexagonal cells with :math:`\sqrt{3}/2` written as 0.866025). Numerical noise then
   produces many spurious constraints. Remedies:

   * Give the lattice vectors and positions with full double precision, or set
     :ref:`SYMMETRIZE_STRUCTURE = 1 <alm_symmetrize_structure>`. *alm* warns when the input deviates from
     the ideal structure of its space group.
   * :ref:`FCSYM_BASIS = Lattice <alm_fcsym_basis>` can be more robust for such structures (it cannot be
     combined with the rotational invariance constraints). The optimal ``L1_ALPHA`` depends on this choice.

   The constraint generation of large models is also expensive by itself. Reduce the cutoff radii or the
   many-body terms (``NBODY``) of the highest orders when it is too slow. *alm* is parallelized with OpenMP;
   check ``Number of OpenMP threads`` at the top of the log.

.. dropdown:: Does alm impose the rotational invariance? The ZA branch of my 2D material is not quadratic.
   :icon: question

   The space-group symmetry and, by default, the translational invariance (acoustic sum rule) are imposed. The
   rotational invariance is imposed with :ref:`ICONST = 12 <alm_iconst>` (or 13, 14) and
   :ref:`ROTAXIS <alm_rotaxis>`, e.g. ``ICONST = 12; ROTAXIS = xyz``. These constraints are formulated in the
   Cartesian frame, need ``FCSYM_BASIS = Cartesian`` (the default), and need a structure symmetric to nearly
   machine precision (see :ref:`SYMMETRIZE_STRUCTURE <alm_symmetrize_structure>`).

   For a 2D material, also make the harmonic supercell large enough in the plane, since the flexural (ZA)
   mode near :math:`\Gamma` is sensitive to the long-range harmonic force constants.

.. dropdown:: alm finds fewer symmetry operations than expected.
   :icon: question

   * Check the structure: the supercell must contain the right number of atoms of each species, and the
     coordinates must have enough digits. ``SYMMETRIZE_STRUCTURE = 1`` symmetrizes a slightly distorted input.
   * The tolerance of the symmetry search is :ref:`TOLERANCE <alm_tolerance>` (default 10\ :sup:`-3`).
   * ``PRINTSYM = 1`` writes the symmetry operations to ``SYMM_INFO``.
   * For magnetic systems, :ref:`MAGMOM <alm_magmom>` lowers the symmetry according to the magnetic order.

.. dropdown:: Can I use snapshots of ab initio molecular dynamics?
   :icon: question

   Yes. The recommended way is to add random displacements to snapshots of the MD trajectory, which reduces
   the correlation between successive snapshots, then calculate their forces with DFT::

      $ displace.py --VASP SPOSCAR0 -md vasprun.xml -e 1001:5000:50 --random --mag 0.04 --prefix disp_random

   The ``-e start:end:interval`` option selects the snapshots. The MD trajectory itself can also be fitted
   after converting it with ``extract.py``, but the resulting force constants are biased toward the MD
   temperature. See :ref:`label_tutorial_bto_ifc`.

   The ``vasprun.xml`` of a VASP machine-learning MD run has a different structure and may not be readable
   by the scripts.


.. _faq_phonons:

Harmonic phonons
----------------

.. dropdown:: The phonon dispersion has imaginary frequencies although the material should be stable.
   :icon: question

   Check the following points.

   1. The structure: relax the primitive cell with the same calculator (DFT functional or force field) that
      you use for the forces, and build the supercell from that relaxed cell. Check that *alm* finds the
      expected space group.
   2. The residual forces and the fitting error (see the questions above), and the displacement magnitude.
   3. The convergence of the DFT forces (cutoff energy, :math:`k` points, SCF threshold).
   4. The supercell size: imaginary modes away from the commensurate :math:`q` points can be an artifact of
      a too small supercell.
   5. Polar materials need the non-analytic correction (``NONANALYTIC`` and ``BORNINFO``).

   Imaginary modes make the lifetimes and thermal conductivity unreliable. They can also be real: the
   structure may be dynamically unstable at the level of theory used, as for high-temperature phases that
   are stabilized by anharmonicity (see :ref:`faq_finite_t`).

.. dropdown:: Why are the phonon dispersion curves discontinuous, or why do I get too many branches?
   :icon: question

   * **Too many (folded) branches.** *anphon* calculates the phonons of the cell in the ``&cell`` field.
     With an HDF5 force-constant file (``PREFIX``.h5), the field can be omitted, and the primitive cell
     stored in the file is used. If a conventional cell or a supercell is given, *anphon* runs on that cell:
     the branches are folded into its smaller Brillouin zone, the :math:`k`-point coordinates refer to its
     reciprocal lattice, and the cost is higher. Note that the lattice vectors in the ``&cell`` field of
     *alm* are those of the supercell.
   * **Jumps at the zone boundaries in polar materials.** With ``NONANALYTIC = 2`` (mixed-space approach),
     the dispersion can jump at zone-boundary points that are not commensurate with the supercell.
     ``NONANALYTIC = 3`` (Ewald method) avoids this (see :ref:`label_tutorial_pbte`).
   * **Discontinuity at** :math:`\Gamma`. In non-cubic polar crystals, the LO-TO splitting depends on the
     direction from which :math:`q` approaches :math:`\Gamma`, so the dispersion is discontinuous there.
     This is physical.
   * **Branches swapped at crossings.** The bands are sorted by frequency by default.
     :ref:`BCONNECT = 1 <anphon_bconnect>` connects them by the similarity of the eigenvectors.

.. dropdown:: Which NONANALYTIC option should I use, and how do I make the BORNINFO file?
   :icon: question

   ``NONANALYTIC = 3`` (Ewald method) is recommended because it separates the short- and long-range parts
   properly. The dielectric tensor and the Born effective charges come from a DFPT calculation, and the file
   format is described in :ref:`label_format_BORNINFO`. For VASP and Quantum ESPRESSO, ``extract.py`` writes
   the file::

      $ extract.py --VASP POSCAR --get born vasprun.xml > BORNINFO

   The Born effective charges must be listed in the order of the atoms of the primitive cell used by
   *anphon*, so the primitive cell of the DFPT calculation must be the same. :ref:`BORNSYM = 1 <anphon_bornsym>`
   symmetrizes the Born charges. The non-analytic correction acts on the harmonic force constants only; it
   is not applied to the cubic and quartic force constants. See the :ref:`PbTe tutorial <label_tutorial_pbte>`
   for a comparison of the options.

.. dropdown:: How should I treat 2D materials and slabs?
   :icon: question

   * In *alm*, :ref:`PERIODIC <alm_periodic>` turns off the periodic boundary condition along the vacuum
     direction when the interacting atoms are searched, e.g. ``PERIODIC = 1 1 0``.
   * The rotational invariance (``ICONST = 12``) and a large in-plane harmonic supercell help to obtain the
     quadratic ZA branch (see above).
   * The thermal conductivity includes the vacuum in the cell volume; see
     :ref:`the question in the thermal conductivity section <faq_kappa_2d>`.

.. dropdown:: How do I get the eigenvectors, animate a mode, or identify the modes at Γ?
   :icon: question

   * :ref:`PRINTEVEC = 1 <anphon_printevec>` writes the eigenvalues and eigenvectors (``PREFIX``.evec.h5 by
     default), in the Cartesian frame.
   * :ref:`ANIME <anphon_anime>` writes animation files of a mode for visualization programs.
   * :ref:`IRREPS <anphon_irreps>` assigns the irreducible representations of the modes at :math:`\Gamma`.
   * :ref:`HESSIAN = 1 <alm_hessian>` in *alm* writes the full force-constant matrix of the supercell
     (``PREFIX``.hessian). It lists the supercell atom and component indices with the force constants only;
     to build the dynamical matrix at arbitrary :math:`q`, combine it with the supercell structure and the masses,
     or simply use *anphon* (or ``anphon.Phonon`` in the :doc:`python_api`).


.. _faq_kappa:

Thermal conductivity
--------------------

.. dropdown:: How do I converge the thermal conductivity?
   :icon: question

   * **q mesh.** Increase the mesh of the ``&kpoint`` field (``KPMODE = 2``) until the thermal conductivity
     converges. The mesh need not be commensurate with the supercell, but it should not break the lattice
     symmetry (e.g. avoid ``10 10 20`` for a cubic cell). Small non-monotonic changes between neighboring
     meshes can occur.
   * **Brillouin-zone integration.** The default is the tetrahedron method (``ISMEAR = -1``). The adaptive
     Gaussian smearing (``ISMEAR = 2``) chooses the widths automatically; with a fixed width
     (``ISMEAR = 0, 1``), check the convergence with respect to ``EPSILON`` as well.
   * **Force constants.** Check the convergence with respect to the cubic cutoff radii and the supercell size.

.. dropdown:: Which force constants are used, and how do I include four-phonon scattering?
   :icon: question

   ``MODE = kappa`` (formerly ``MODE = RTA``) needs the harmonic and cubic force constants; higher orders in
   the file are not used for the three-phonon calculation. Four-phonon scattering is added with
   :ref:`INCLUDE_4PH = 1 <anphon_include_4ph>`, which needs quartic force constants. Because it is much more
   expensive, the four-phonon rates can be evaluated on a coarser mesh (:ref:`KMESH_COARSE <anphon_kmesh_coarse>`)
   and interpolated. See the :ref:`BAs tutorial <label_tutorial_bas_4ph>`.

.. dropdown:: Can the BTE be solved beyond the relaxation-time approximation?
   :icon: question

   Yes. :ref:`SOLVER <anphon_solver>` in the ``&kappa`` field selects the solver: ``RTA`` (default),
   ``IBTE`` (iterative), ``VBTE`` (variational), and ``DBTE`` (direct, for small meshes). The results of the
   iterative solvers are written to ``PREFIX``.kl_iter. See :ref:`kappa_beyond_rta` for the theory.

.. dropdown:: Is isotope or boundary scattering included?
   :icon: question

   * Isotope scattering is included by default (:ref:`ISOTOPE = 1 <anphon_isotope>`) with the natural
     abundances of the elements. Set ``ISOTOPE = 0`` for isotopically pure crystals, or give ``ISOFACT``
     for enriched samples.
   * Boundary scattering is added with :ref:`LEN_BOUNDARY <anphon_len_boundary>`, or afterwards with
     ``analyzer.py --calc kappa --size L --h5 PREFIX.kappa.h5`` (``L`` in nm).

.. dropdown:: How do I get the phonon lifetimes, mean free paths, or the cumulative thermal conductivity?
   :icon: question

   ``analyzer.py`` reads ``PREFIX``.kappa.h5::

      $ analyzer.py --calc tau --temp 300 --h5 PREFIX.kappa.h5
      $ analyzer.py --calc cumulative --temp 300 --h5 PREFIX.kappa.h5

   The first command prints the lifetime (ps), group velocity (m/s), and mean free path (nm) of each mode.
   The layout of ``PREFIX``.kappa.h5 is described in :ref:`label_hdf5_kappa`, and ``KAPPA_SPEC = 1`` writes
   the spectrum of the thermal conductivity. For the linewidth or self-energy of selected modes at arbitrary
   :math:`q` points or along a band path, use ``MODE = selfenergy`` (see the
   :ref:`&selfenergy field <anphon_selfenergy_field>`).

.. _faq_kappa_2d:

.. dropdown:: The thermal conductivity of my 2D material is much smaller than reported values.
   :icon: question

   The thermal conductivity is divided by the volume of the cell, which includes the vacuum. Multiply it by
   :math:`c/d`, where :math:`c` is the length of the cell perpendicular to the sheet and :math:`d` the
   thickness of the sheet, and use the same :math:`d` when comparing with other works.
   The thickness :math:`d` is a convention that must be stated together with the result (e.g., 3.35 Å, the interlayer
   distance of graphite, is commonly used for graphene).
   See `discussion #32 <https://github.com/alamode-code/alamode/discussions/32>`_.


.. _faq_finite_t:

Finite temperature (SCPH and QHA)
---------------------------------

.. dropdown:: The SCPH iteration does not converge.
   :icon: question

   * Use a smaller :ref:`MIXALPHA <anphon_mixalpha>`.
   * Use a smaller temperature step ``DT``: with ``WARMSTART = 1`` (default), each temperature starts from the
     solution of the previous one, so a closer starting point helps. ``LOWER_TEMP`` sets the direction of the
     temperature sweep.
   * Check the quartic force constants: inaccurate ones (too small cutoff radii, too few training data) are a
     common cause. Check how the results change when you add training data.

.. dropdown:: Imaginary modes remain (or appear) after the SCPH calculation.
   :icon: question

   * The SCPH equation is solved on the :math:`q` points of :ref:`KMESH_INTERPOLATE <anphon_kmesh_interpolate>`
     only; the frequencies elsewhere are interpolated. Imaginary modes are therefore detected only on this
     mesh, and with ``KMESH_INTERPOLATE = 1 1 1`` only at :math:`\Gamma`. Imaginary pockets can appear at the
     interpolated points.
   * Small imaginary pockets near :math:`\Gamma` often indicate that the quartic force constants are not
     accurate enough. Increase the quartic cutoff radius and the number of training data.
   * Use the same ``NONANALYTIC`` option in the harmonic and SCPH calculations.
   * The phonons can be genuinely unstable at low temperatures, e.g. below a phase transition.

.. dropdown:: How do I choose KMESH_INTERPOLATE and KMESH_SCPH?
   :icon: question

   The :math:`q` points of ``KMESH_INTERPOLATE`` should be commensurate with the supercell of the harmonic
   force constants: for a 2x2x2 supercell, use ``KMESH_INTERPOLATE = 2 2 2`` (or ``1 1 1``, which is less
   accurate). anphon can evaluate the harmonic phonons on any mesh, but the anharmonic renormalization is
   interpolated from the ``KMESH_INTERPOLATE`` points, so converge this mesh and the supercell size
   (the range of the force constants) separately. ``KMESH_SCPH`` can be denser, but each of its values must be a multiple of the corresponding value of
   ``KMESH_INTERPOLATE``. See the :ref:`SrTiO3 tutorial <label_tutorial_sto_scph>`.

.. dropdown:: How do I use the SCPH (or QHA) phonons for the thermal conductivity, PDOS, or eigenvectors?
   :icon: question

   The renormalized harmonic force constants at each temperature are stored in ``PREFIX``.scph.h5
   (``PREFIX``.qha.h5). In a later *anphon* run (``MODE = phonons`` or ``kappa``), keep ``FCSFILE`` pointing
   to the original force constants and give

   .. code-block:: none

      DFC2FILE = PREFIX.scph.h5
      FC2_TEMPERATURE = 300

   in the ``&general`` field. All harmonic properties (bands, DOS, PDOS, eigenvectors) are then those of the
   SCPH phonons at 300 K. For ``MODE = kappa``, also set ``TMIN = TMAX = FC2_TEMPERATURE`` and repeat the run
   for each temperature. The harmonic force constants may come from a larger supercell than the SCPH run.
   If the structure was relaxed (``RELAX_STR > 0``), add ``RELAXED_STRUCTURE = 1``.
   Use the same ``NONANALYTIC`` (2 and 3 are interchangeable) and ``BORNINFO`` as in the SCPH run;
   other combinations stop with an error or print a warning. The SCPH frequencies are reproduced exactly
   at the :math:`q` points of ``KMESH_INTERPOLATE``; elsewhere they are interpolated.
   This replaces the former ``dfc2`` workflow. See :ref:`FC2_TEMPERATURE <anphon_fc2_temperature>` and
   :ref:`DFC2FILE <anphon_dfc2file>`.

.. dropdown:: My SCPH run finishes immediately and ignores the new settings.
   :icon: question

   When ``PREFIX``.scph.h5 exists in the working directory, ``RESTART_SCPH`` defaults to 1 and anphon restarts
   from the stored solution instead of solving the SCPH equation again. If ``KMESH_INTERPOLATE``,
   ``KMESH_SCPH``, or ``NONANALYTIC`` differ from the stored run, it stops with an error instead.
   Delete (or rename) the file, use a new ``PREFIX``, or set ``RESTART_SCPH = 0`` to recalculate.
   (``BUBBLE = 4`` never restarts automatically.)

.. dropdown:: What do the BUBBLE options do?
   :icon: question

   :ref:`BUBBLE <anphon_bubble>` adds the bubble (three-phonon) self-energy on top of the SCPH solution. The
   options 1, 2, and 3 differ in the frequency at which the self-energy is evaluated (zero, the SCPH
   frequency, or self-consistently the corrected frequency). The correction can lower the frequencies strongly
   and even make a stable mode imaginary; it is a frequency correction, not a stability test. ``BUBBLE = 4``
   evaluates the curvature of the free energy at :math:`\Gamma` of the chosen cell for the stability question;
   it runs along a structural optimization and needs ``RELAX_STR != 0``, ``KMESH_SCPH = KMESH_INTERPOLATE``,
   and ``NONANALYTIC = 0``, and it cannot be combined with ``RESTART_SCPH = 1``.
   For options 1--3, the :math:`q` mesh of the bubble sum is set by :ref:`KMESH_BUBBLE <anphon_kmesh_bubble>`.

.. dropdown:: The thermal conductivity calculated with the SCPH phonons increases with temperature. Is it wrong?
   :icon: question

   Not necessarily. When the phonons harden with temperature, the larger frequencies and group velocities
   can outweigh the increase of the scattering rates. If the SCPH frequencies are lower than the harmonic
   ones, the opposite happens. See
   `discussion #319 <https://github.com/alamode-code/alamode/discussions/319>`_.

.. dropdown:: How do I optimize the crystal structure at finite temperatures?
   :icon: question

   :ref:`RELAX_STR <anphon_relax_str>` (SCPH) and the corresponding tag of the ``&qha`` field relax the
   atomic positions, or the positions and the cell, at each temperature. Relaxing the cell needs the
   coupling between strain and force constants (see :ref:`label_strain_tools` and
   :ref:`STRAIN_COUPLING <anphon_strain_coupling>`). The relaxed structures are stored in ``PREFIX``.scph.h5
   / ``PREFIX``.qha.h5. See the tutorials :ref:`label_tutorial_bto_scph_relax` and
   :ref:`label_tutorial_zno_qha_relax`.


.. _faq_interfaces:

Interfaces, file conversion, and performance
--------------------------------------------

.. dropdown:: Can I use a DFT code that displace.py and extract.py do not support?
   :icon: question

   Yes. ``displace.py`` and ``extract.py`` support VASP, Quantum ESPRESSO, OpenMX, xTAPP, and LAMMPS, but
   *alm* only needs a ``DFSET`` file with the displacements and forces of each structure. :doc:`other_codes`
   shows how to make it with ASE calculators and what to look out for with codes such as WIEN2k and CRYSTAL.

.. dropdown:: How do I convert force constants to or from phonopy, phono3py, or ShengBTE?
   :icon: question

   Use ``tools/fcconvert.py``, which converts harmonic and cubic force constants between ALAMODE, phonopy/phono3py,
   and ShengBTE, and quartic force constants between ALAMODE and FourPhonon, in both directions; see :ref:`label_fc_conversion`. *alm* can also write the cubic and quartic force
   constants for ShengBTE and FourPhonon directly with :ref:`FC3_SHENGBTE <alm_fc3_shengbte>` and
   :ref:`FC4_SHENGBTE <alm_fc4_shengbte>`.

   When you combine data of different codes, the atoms of the supercells must correspond one to one: atom
   :math:`i` of one code must be atom :math:`i` of the other. A mismatch gives imaginary frequencies or a
   wrong thermal conductivity. See
   `discussion #301 <https://github.com/alamode-code/alamode/discussions/301>`_.

.. dropdown:: Can I use harmonic force constants from DFPT (Quantum ESPRESSO)?
   :icon: question

   Yes. ``qe2alm`` converts a Quantum ESPRESSO force-constant file (output of ``q2r.x``) to the ALAMODE XML
   format, which *anphon* reads. It asks for the input and output file names interactively. Conversely,
   ``scph_to_qefc.py`` writes the SCPH-renormalized harmonic force constants in the Quantum ESPRESSO format.
   A DFPT calculation on an :math:`N\times N\times N` :math:`q` mesh and a finite-displacement calculation in
   the :math:`N\times N\times N` supercell give almost the same harmonic force constants when the same code
   and parameters are used and all harmonic interactions are included (cutoff ``None``).

.. dropdown:: How do I read the force constants in my own scripts?
   :icon: question

   ``PREFIX``.h5 is an HDF5 file; its layout is described in :ref:`label_hdf5_fcs`. ``fcsinfo.py`` prints a
   summary of its contents. The Python package ``alm`` (:ref:`label_python_api`) returns the force constants as
   NumPy arrays.

.. dropdown:: alm stops with std::bad_alloc.
   :icon: question

   The default dense solver stores the sensing matrix, whose size is about
   (number of structures) :math:`\times` (3 :math:`\times` number of atoms) :math:`\times` (number of
   parameters) :math:`\times` 8 bytes. Options:

   * ``LMODEL = ols``: :ref:`SPARSE = 1 <alm_sparse>` for data in which only a few atoms are displaced
     (the zero displacements must be exactly zero in ``DFSET``), or :ref:`USE_CHOLESKY = 1 <alm_use_cholesky>`,
     which builds :math:`A^{\mathsf T}A` in chunks (:ref:`CHUNKSIZE <alm_chunksize>`). ``SPARSE = 1`` has no
     effect for the elastic net and adaptive LASSO.
   * ``LMODEL = enet`` or ``adaptive-lasso``: :ref:`L1_SOLVER = fista <alm_l1_solver>` avoids the dense Gram
     matrix.
   * Reduce the model: smaller cutoff radii, fewer many-body terms (``NBODY``), or fewer structures.

.. dropdown:: How are alm and anphon parallelized?
   :icon: question

   * *alm* uses OpenMP only (``OMP_NUM_THREADS``); it has no MPI parallelization.
   * *anphon* uses MPI and OpenMP, which can be combined. For SCPH with many atoms in the primitive cell and
     few :math:`k` points, :ref:`IALGO = 1 <anphon_ialgo>` distributes the phonon branches instead of the
     :math:`k` points. :ref:`V4_REAL_SPACE = 1 <anphon_v4_real_space>` avoids storing the quartic coupling
     tensor, whose memory grows as the fourth power of the number of modes; it runs on one MPI process with
     many OpenMP threads.
