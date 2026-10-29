ANPHON: Theoretical background
==============================

.. |umulaut_u|   unicode:: U+00FC
.. |umulaut_o|   unicode:: U+00F6

Dynamical matrix
----------------

The dynamical matrix is given by

.. math::
    :label: dymat 

    D_{\mu\nu}(\kappa\kappa^{\prime};\boldsymbol{q}) = \frac{1}{\sqrt{M_{\kappa}M_{\kappa^{\prime}}}}
    \sum_{\ell^{\prime}}\Phi_{\mu\nu}(\ell\kappa;\ell^{\prime}\kappa^{\prime})\exp{\left[i\boldsymbol{q}\cdot(\boldsymbol{r}(\ell^{\prime})-\boldsymbol{r}(\ell))\right]},

where :math:`M_{\kappa}` is the atomic mass of atom :math:`\kappa`.
By diagonalizing the dynamical matrix, one can obtain :math:`m` (:math:`=3N_{\kappa}`) eigenvalues :math:`\omega_{\boldsymbol{q}j}^{2}`  (:math:`j = 1, 2, \dots, m`) and corresponding eigenvectors :math:`\boldsymbol{e}_{\boldsymbol{q}j}` for each :math:`\boldsymbol{q}` point.
Here, :math:`\boldsymbol{e}_{\boldsymbol{q}j}` is a column vector consisting of atomic polarization :math:`e_{\mu}(\kappa;\boldsymbol{q}j)`.
Let :math:`D(\boldsymbol{q})` denote a matrix form of equation :eq:`dymat`, the eigenvalues may be written as

.. math::
    :label: omega2

    \omega_{\boldsymbol{q}j}^{2} = (\boldsymbol{e}_{\boldsymbol{q}j}^{*})^{\mathrm{T}} D(\boldsymbol{q})\boldsymbol{e}_{\boldsymbol{q}j}.

Next, we introduce :math:`m\times m` matrices :math:`\Lambda` and :math:`W` which are defined as 
:math:`\Lambda(\boldsymbol{q}) = \mathrm{diag} (\omega_{\boldsymbol{q}1}^{2},\dots,\omega_{\boldsymbol{q}m}^{2})` and 
:math:`W(\boldsymbol{q}) = (\boldsymbol{e}_{\boldsymbol{q}1},\dots,\boldsymbol{e}_{\boldsymbol{q}m})`, respectively. 
Then, equation :eq:`omega2` can be denoted as 

.. math::
    
    \Lambda(\boldsymbol{q}) = W^{\dagger}(\boldsymbol{q})D(\boldsymbol{q})W(\boldsymbol{q}).


Non-analytic correction
-----------------------

When one needs to capture the LO\-TO splitting near the zone-center, it is necessary to account for the long\-range correction to the dynamical matrix. 
In the **anphon** code, three different approaches are implemented: the Parlinski's method [Parlinski1998]_ (``NONANALYTIC = 1``), the mixed-space approach [Wang2010]_  (``NONANALYTIC = 2``), and the Ewald summation method (``NONANALYTIC = 3``). 

The first two methods are rather simple and add the non-analytic part of the dynamical matrix with different weighting factors.
The non-analytic part of the dynamical matrix is defined as

.. math::

    D_{\mu\nu}^{\mathrm{NA}}(\kappa\kappa^{\prime};\boldsymbol{q}) = \frac{1}{\sqrt{M_{\kappa}M_{\kappa^{\prime}}}}
    \frac{4\pi e^{2}}{V} \frac{(Z_{\kappa}^{*}\boldsymbol{q})_{\mu}(Z_{\kappa^{\prime}}^{*}\boldsymbol{q})_{\nu}}{\boldsymbol{q}\cdot\epsilon^{\infty}\boldsymbol{q}},

where :math:`V` is the volume of the primitive cell, :math:`Z_{\kappa}^{*}` is the Born effective charge tensor of atom :math:`\kappa`, 
and :math:`\epsilon^{\infty}` is the electronic dielectric permittivity tensor, respectively.

Parlinski's approach (``NONANALYTIC = 1``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the Parlinski's approach [Parlinski1998]_ , the total dynamical matrix is given by

.. math::

    D(\boldsymbol{q}) + D^{\textrm{NA}}(\boldsymbol{q})\exp{(-q^{2}/\sigma^{2})},

where :math:`\sigma` is a damping factor (``NA_SIGMA``). 
:math:`\sigma` must be chosen carefully so that the non-analytic contribution becomes negligible at Brillouin zone (BZ) boundaries. 
Too large :math:`\sigma` can result in a discontinuity of phonon dispersion at the BZ boundaries, while too small :math:`\sigma` can yield a large oscillation of the LO phonon frequencies as going away from :math:`\boldsymbol{q}=0`. So, an appropriate :math:`\sigma` value needs to be chosen in an *ad hoc* way.


Mixed-space approach (``NONANALYTIC = 2``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the mixed-space approach [Wang2010]_ , the total dynamical matrix is given by

.. math::

    D(\boldsymbol{q}) + D^{\textrm{NA}}(\boldsymbol{q})\frac{1}{N}\sum_{\ell^{\prime}}\exp{\left[i\boldsymbol{q}\cdot(\boldsymbol{r}(\ell^{\prime})-\boldsymbol{r}(\ell))\right]}.

The second term vanishes at the :math:`\boldsymbol{q}` points commensurate with the supercell except for the :math:`\Gamma` point (:math:`\boldsymbol{q} = 0`).
No adjustable parameters exist in this approach. In many cases, the mixed-space approach gives reasonable dispersion curves. 
However, it can still yield a discontinuity in the phonon dispersion at BZ boundaries when the boundary point is incommensurate with the supercell size.


Ewald method (``NONANALYTIC = 3``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Ewald summation technique [Gonze1997]_ is more accurate and therefore recommended in general.
In this method, we first decompose the harmonic force constants computed in a supercell into the long\-range dipole-dipole (DD) components and the short\-range ones as follows

.. math:: 

    \Phi_{\mu\nu}^{\mathrm{SR}}(\ell\kappa;\ell^{\prime}\kappa^{\prime})= \Phi_{\mu\nu}(\ell\kappa;\ell^{\prime}\kappa^{\prime}) - \Phi_{\mu\nu}^{\mathrm{DD}}(\ell\kappa;\ell^{\prime}\kappa^{\prime}).


After the short\-range components are obtained, they are used in equation :eq:`dymat` to obtain the short\-range dynamical matrix :math:`D_{\mu\nu}^{\mathrm{SR}}(\kappa\kappa^{\prime};\boldsymbol{q})`. Then, the total dynamical matrix is constructed as

.. math:: 

    D_{\mu\nu}(\kappa\kappa^{\prime};\boldsymbol{q}) = D^{\mathrm{SR}}_{\mu\nu}(\kappa\kappa^{\prime};\boldsymbol{q}) +  D^{\mathrm{DD}}_{\mu\nu}(\kappa\kappa^{\prime};\boldsymbol{q}).

The long\-range terms, :math:`\Phi_{\mu\nu}^{\mathrm{DD}}(\ell\kappa;\ell^{\prime}\kappa^{\prime})` and :math:`D^{\mathrm{DD}}_{\mu\nu}(\kappa\kappa^{\prime};\boldsymbol{q})`, can be calculated from :math:`\epsilon^{\infty}` and :math:`Z_{\kappa}^{*}` using the Ewald summation technique.

To include the non-analytic correction with ``NONANALYTIC > 0``, one also needs to give ``BORNINFO``.


Irreducible representations and IR/Raman activity (``IRREPS = 1``)
------------------------------------------------------------------

At the zone center, every space-group operation :math:`\{S|\boldsymbol{v}\}` acts on the
mass-weighted polarization vectors through the :math:`3N_{\kappa}\times 3N_{\kappa}` matrix

.. math::
   T_{\mu\nu}(S;\kappa^{\prime}\kappa) = S_{\mu\nu}\,\delta_{\kappa^{\prime},S(\kappa)},

where :math:`S(\kappa)` is the atom onto which :math:`\kappa` is mapped; the fractional
translation :math:`\boldsymbol{v}` contributes no phase at :math:`\boldsymbol{q}=0`, so the
representation is an ordinary representation of the crystal point group.
For each set :math:`\lambda` of degenerate modes, the character

.. math::
   \chi_{\lambda}(S) = \sum_{j\in\lambda}
   \boldsymbol{e}_{j}^{\dagger}\, T(S)\, \boldsymbol{e}_{j}

is a class function, and the multiplicity of the irreducible representation
:math:`\mu` in :math:`\lambda` follows from the orthogonality relation

.. math::
   n_{\mu}(\lambda) = \frac{1}{g\,c_{\mu}} \sum_{S} \chi_{\lambda}(S)\,\chi_{\mu}(S)^{*},

where :math:`g` is the order of the point group and :math:`c_{\mu}` accounts for the
time-reversal pairing of complex-conjugate one-dimensional representations
(:math:`c_{\mu}=2` for the paired :math:`E`-type representations of the groups
:math:`C_{3}`, :math:`C_{4}`, :math:`C_{6}`, :math:`S_{4}`, :math:`S_{6}`, :math:`C_{3h}`,
:math:`C_{4h}`, :math:`C_{6h}`, :math:`T`, :math:`T_{h}`, and :math:`c_{\mu}=1` otherwise).
A multiplet is infrared active when its representation appears in the vector
representation :math:`\Gamma_{V}` (characters :math:`\operatorname{tr}S`), and Raman
active when it appears in the symmetric square :math:`[\Gamma_{V}\otimes\Gamma_{V}]_{s}`
(characters :math:`(\operatorname{tr}^{2}S + \operatorname{tr}S^{2})/2`).
The analysis uses the analytic dynamical matrix, i.e. the labels refer to the
:math:`\boldsymbol{q}\rightarrow 0` limit without the direction-dependent non-analytic
correction.

When ``BORNINFO`` is given, the infrared oscillator-strength tensor of a multiplet is
also computed from the mode effective charges
:math:`Z^{*}_{\nu,\alpha} = \sum_{\kappa\beta} Z^{*}_{\kappa,\alpha\beta}\,
e_{\beta}(\kappa;\nu)/\sqrt{M_{\kappa}}` as

.. math::
   S_{\alpha\beta}(\lambda) = \sum_{\nu\in\lambda} Z^{*}_{\nu,\alpha} Z^{*}_{\nu,\beta},

which is invariant under rotations within the degenerate subspace.

.. _group_velocity:

Group velocity
--------------

The group velocity of phonon mode :math:`\boldsymbol{q}j` is given by 

.. math::
    
    \boldsymbol{v}_{\boldsymbol{q}j} = \frac{\partial \omega_{\boldsymbol{q}j}}{\partial \boldsymbol{q}}.

*anphon* evaluates it from the **velocity matrix** (default) or by a **finite difference**.

Velocity matrix (default)
~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`\tilde{D}(\boldsymbol{q})` denote the dynamical matrix whose phase factor carries the full interatomic
separation rather than the lattice vector alone,

.. math::
    :label: dymat_tilde

    \tilde{D}_{\mu\nu}(\kappa\kappa^{\prime};\boldsymbol{q}) = \frac{1}{\sqrt{M_{\kappa}M_{\kappa^{\prime}}}}
    \sum_{\ell^{\prime}}\Phi_{\mu\nu}(\ell\kappa;\ell^{\prime}\kappa^{\prime})
    \exp{\left[i\boldsymbol{q}\cdot(\boldsymbol{r}(\ell^{\prime}\kappa^{\prime})-\boldsymbol{r}(\ell\kappa))\right]}.

It is related to :eq:`dymat` by the unitary transformation
:math:`\tilde{D}(\boldsymbol{q}) = U^{\dagger}(\boldsymbol{q})D(\boldsymbol{q})U(\boldsymbol{q})` with
:math:`U_{\kappa\kappa^{\prime}}(\boldsymbol{q})=\delta_{\kappa\kappa^{\prime}}\,e^{i\boldsymbol{q}\cdot\boldsymbol{r}(\kappa)}`,
so it has the same eigenvalues :math:`\omega_{\boldsymbol{q}j}^{2}`, and its eigenvectors are
:math:`\tilde{\boldsymbol{e}}_{\boldsymbol{q}j}=U^{\dagger}(\boldsymbol{q})\boldsymbol{e}_{\boldsymbol{q}j}`.
The band off-diagonal generalization of the group velocity [Allen1993]_ [Simoncelli2022]_ is

.. math::
    :label: velmat

    v_{\boldsymbol{q}jj'}^{\mu} = \frac{1}{2\sqrt{\omega_{\boldsymbol{q}j}\omega_{\boldsymbol{q}j'}}}
    (\tilde{\boldsymbol{e}}_{\boldsymbol{q}j}^{*})^{\mathrm{T}}
    \frac{\partial \tilde{D}(\boldsymbol{q})}{\partial q_{\mu}} \tilde{\boldsymbol{e}}_{\boldsymbol{q}j'},

whose diagonal :math:`v_{\boldsymbol{q}jj}^{\mu}` is the group velocity. It is :math:`\tilde{D}`, not :math:`D`, that must
be differentiated: the two derivatives differ by :math:`i\,[\,D,\,\mathrm{diag}(r_{\mu}(\kappa))\,]`, whose matrix elements
between eigenvectors carry the factor :math:`\omega_{\boldsymbol{q}j}^{2}-\omega_{\boldsymbol{q}j'}^{2}`. The diagonal and any
element inside a degenerate multiplet are therefore the same for both, but genuinely off-diagonal elements --- those entering the
:ref:`coherent term <kappa_coherent>` --- are not (the *displacement-aware* convention of Ref. [Simoncelli2022]_).

For polar systems the non-analytic correction (``NONANALYTIC = 1, 2, 3``) contributes to
:math:`\partial\tilde{D}/\partial\boldsymbol{q}`. It is obtained by a central difference of the assembled non-analytic
matrix together with the derivative of the sublattice phase factor relating :math:`D` and :math:`\tilde{D}`; omitting the
latter corrupts even the diagonal velocities. At :math:`\Gamma` the non-analytic velocity is not defined without a
directional convention and is set to zero. On a band path the derivative uses the same direction vector as the eigenproblem
for that segment.

At exactly :math:`\boldsymbol{q}=0` the three acoustic modes are uniform translations whose frequencies vanish by the acoustic
sum rule, and :math:`v_{\boldsymbol{q}jj'}^{\mu}` with :math:`j` or :math:`j'` among them is :math:`0/0`. Velocity-matrix elements
are set to zero wherever a frequency is below :math:`10^{-8}` Ry, but a small sum-rule residual can leave these modes just above it,
and their elements are then amplified without bound (e.g. in polar crystals without inversion). The translational modes are
therefore identified at :math:`\Gamma` as the three eigenvectors with the largest weight in the span of the mass-weighted
uniform translations :math:`t_{\alpha}(\kappa i)\propto\sqrt{m_{\kappa}}\,\delta_{i\alpha}`, independently of their frequencies,
and their rows and columns of the velocity matrix are set to zero, including the coherent elements with optical modes; their
finite-difference velocities (zero by time reversal up to roundoff) are set to zero as well. This applies to the transport terms
and to ``PRINTVEL``. It is a convention for the single point :math:`\boldsymbol{q}=0`, not a sound velocity: sound velocities follow
from :math:`\boldsymbol{q}\neq0`. The adaptive smearing widths are left unchanged.

Inside a degenerate multiplet :math:`\mathcal{B}` the eigenvectors are fixed only up to a unitary rotation, so the individual
:math:`v_{\boldsymbol{q}jj}^{\mu}` are not defined; only block traces such as
:math:`\sum_{j,j'\in\mathcal{B}} v_{\boldsymbol{q}jj'}^{\mu}v_{\boldsymbol{q}j'j}^{\nu}` are. How the transport terms use this is described
in :ref:`the Peierls-term section <kappa_peierls>`.

On a uniform mesh the velocity operator is averaged over the little group :math:`G_{\boldsymbol{q}}` of
:math:`\boldsymbol{q}` before it is projected onto the eigenvectors: over the operations :math:`R` of the space group used for
the k-point reduction with :math:`R\boldsymbol{q}=\boldsymbol{q}+\boldsymbol{G}`, and, combined with time reversal, those with
:math:`R\boldsymbol{q}=-\boldsymbol{q}+\boldsymbol{G}`. With :math:`M^{\mu}=\partial\tilde{D}/\partial q_{\mu}` in the atomic
Cartesian basis and :math:`T(R)` the transformation of the displacement field (atom permutation, rotation and cell phase),

.. math::

    M^{\mu} \leftarrow \frac{1}{|G_{\boldsymbol{q}}|}\sum_{R\in G_{\boldsymbol{q}}}\sum_{\nu}R_{\mu\nu}\,
    T(R)\,M^{\nu}\,T(R)^{\dagger},

where an antiunitary element enters as :math:`-[T(R)M^{\nu}T(R)^{\dagger}]^{*}`. The :math:`T(R)` form a representation of
:math:`G_{\boldsymbol{q}}` up to phases, which cancel, so this is an exact projector that does not involve the eigenvectors or
the degeneracies, and it leaves an operator that already has the symmetry unchanged. It matters with ``NONANALYTIC = 1, 2``: their
non-analytic term depends on the direction of :math:`\boldsymbol{q}` itself, so the dynamical matrix is not periodic, and on the
zone boundary the velocities would depend on which of the equally short images :math:`\boldsymbol{q}+\boldsymbol{G}` the mesh
keeps and break the crystal symmetry (e.g. a nonzero :math:`\kappa_{xy}` in a cubic crystal). The eigenvectors at such a point
are those of the one image used and are not symmetrized.

Two limitations remain. First, the average acts on each mesh point separately and does not make the velocities of different
members of a star consistent with each other. Where all equally short images of :math:`\boldsymbol{q}` are related by
:math:`G_{\boldsymbol{q}}` (simple cubic and tetragonal lattices) this follows; where they are not, as on the edges where two
hexagonal faces of the fcc Brillouin zone meet, an anisotropy of order :math:`10^{-6}` of :math:`\kappa` can remain with
``NONANALYTIC = 1, 2`` (e.g. zincblende BAs on a :math:`20\times20\times20` mesh). Second, the statement that the conductivity
uses only block traces is exact only for exactly degenerate blocks: blocks are formed with a tolerance of
:math:`10^{-6}` cm\ :sup:`-1`, the heat capacity is evaluated with each branch's own frequency, and the coherent term uses the
individual frequencies of the pair.

Finite difference
~~~~~~~~~~~~~~~~~

The finite-difference group velocity is obtained from a central difference,

.. math::

    \boldsymbol{v}_{\boldsymbol{q}j} \approx \frac{\omega_{\boldsymbol{q}+\Delta\boldsymbol{q}j} - \omega_{\boldsymbol{q}-\Delta\boldsymbol{q}j}}{2\Delta\boldsymbol{q}},

where :math:`j` is the index in the *sorted* eigenvalue list at each shifted point. Wherever branches cross or are degenerate the sorted
index exchanges character between :math:`\boldsymbol{q}\pm\Delta\boldsymbol{q}` and the quotient connects two different branches.
Finite-difference velocities are used by the adaptive smearing widths (``ISMEAR = 2``) and by the iterative Boltzmann solvers, which
have not been reformulated. At a mesh point with a nontrivial little group they are symmetrized as well: inside a
degenerate block the central difference of sorted eigenvalues pairs the branches and is not a vector, so each block is first
replaced by its mean, which approximates :math:`\mathrm{Tr}(P_{\mathcal{B}}v^{\mu}P_{\mathcal{B}})/d`, and the result is projected
onto the invariant vectors, :math:`\boldsymbol{v}\leftarrow|G_{\boldsymbol{q}}|^{-1}\sum_{R}(\pm R)\boldsymbol{v}` (minus sign
for the time-reversal partners). The block mean is a trace approximation and is used only where the velocity enters
transport as a vector (the iterative and direct Boltzmann solvers and the velocities written to the result files). The
adaptive smearing widths keep the unaveraged differences, because they are built from :math:`|\boldsymbol{v}_{j}\pm\boldsymbol{v}_{j'}|`
per branch, and two opposite slopes in a block give a nonzero width although their mean vanishes; the 3ph linewidths are
therefore unchanged. At a mesh point with a trivial little group nothing is done, so an accidental degeneracy there keeps the
pairing of the sorted branches.

If one needs to save the group velocities, please turn on the ``PRINTVEL``-tag; the printed values follow the same formulation as the
conductivity in that run. At a degeneracy they are one admissible basis choice, not a unique value.


Mode effective charge
---------------------

The mode effective charge is defined as 

.. math::
    
    Z^{*}_{j,\alpha} = \sum_{\kappa\beta}Z^{*}_{\kappa,\alpha\beta}\frac{e_{\beta}(\kappa;\boldsymbol{0}j)}{\sqrt{M_{\kappa}}}

where :math:`Z^{*}_{\kappa,\alpha\beta}` is the atomic Born effective charge. To compute the mode effective charges of the zone-center modes, please set ``ZMODE = 1`` in the &analysis field and supply ``BORNINFO``.


Thermodynamic functions
-----------------------

The specific heat at constant volume :math:`C_{\mathrm{v}}`, the internal energy :math:`U`, 
the vibrational entropy :math:`S`, and the Helmholtz free energy :math:`F` of independent harmonic oscillators are
given as follows:

.. math::
    :nowrap:
    
    \begin{align*}
     U &= \frac{1}{N_{q}}\sum_{\boldsymbol{q},j} \hbar\omega_{\boldsymbol{q}j} \left[\frac{1}{e^{\hbar\omega_{\boldsymbol{q}j}/kT} - 1} + \frac{1}{2}\right], \\
     C_{\mathrm{v}} &= \frac{k}{N_{q}}\sum_{\boldsymbol{q},j} \left(\frac{\hbar\omega_{\boldsymbol{q}j}}{2kT}\right)^{2} \mathrm{cosech}^{2}\left(\frac{\hbar\omega_{\boldsymbol{q}j}}{2kT}\right),\\
     S &= \frac{k}{N_{q}}\sum_{\boldsymbol{q},j} \left[\frac{\hbar\omega_{\boldsymbol{q}j}}{kT} \frac{1}{e^{\hbar\omega_{\boldsymbol{q}j}/kT} - 1} 
        - \log{\left( 1 - e^{-\hbar\omega_{\boldsymbol{q}j}/kT}\right)}\right], \\
     F &= \frac{1}{N_{q}}\sum_{\boldsymbol{q},j}\left[ \frac{\hbar\omega_{\boldsymbol{q}j}}{2} + kT\log{\left( 1 - e^{-\hbar\omega_{\boldsymbol{q}j}/kT}\right)} \right].
    \end{align*}

Here, :math:`k` is the Boltzmann constant. These quantities are saved in the ``PREFIX``.thermo file.

When the self-consistent phonon mode (``MODE = SCPH``) is selected, the anharmonic free-energy 
defined by the following equation will be calculated and saved in the ``PREFIX``.scph_thermo file:

.. math::
    :nowrap:
    
    \begin{align*}
     F^{\mathrm{SCP}} &= \frac{1}{N_{q}}\sum_{\boldsymbol{q},j}\left[ \frac{\hbar\Omega_{\boldsymbol{q}j}}{2} + kT\log{\left( 1 - e^{-\hbar\Omega_{\boldsymbol{q}j}/kT}\right)} \right] \\
     & - \frac{1}{4N_{q}}\sum_{\boldsymbol{q},j}\left[ \Omega_{\boldsymbol{q}j}^{2} - (C_{\boldsymbol{q}}^{\dagger}\Lambda_{\boldsymbol{q}}^{(\mathrm{HA})}C_{\boldsymbol{q}})_{jj} \right]
     \times \frac{\hbar [1 + 2n_{\boldsymbol{q}j} ]}{2\Omega_{\boldsymbol{q}j}}.
    \end{align*}

Details of the derivation of the above expression can be found in Ref. [Oba2019]_.

.. _formalism_fe_bubble:

Bubble correction to the free energy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``FE_BUBBLE = 1``, the lowest-order contribution of the cubic anharmonicity to the free energy, which
comes from the bubble diagram with two cubic vertices, is also computed:

.. math::
    :label: fe_bubble

    F^{\mathrm{B}} &= -\frac{1}{6\hbar N_{q}}\sum_{\boldsymbol{q},\boldsymbol{q}_{1},\boldsymbol{q}_{2}}\sum_{j,j_{1},j_{2}}
    |V^{(3)}_{\boldsymbol{q}j,\boldsymbol{q}_{1}j_{1},\boldsymbol{q}_{2}j_{2}}|^{2} \notag \\
    & \times \left[ \frac{(1+n)(1+n_{1})(1+n_{2}) - n n_{1} n_{2}}{\omega + \omega_{1} + \omega_{2}}
    + 3\,\frac{n(1+n_{1})(1+n_{2}) - (1+n) n_{1} n_{2}}{\omega_{1} + \omega_{2} - \omega} \right],

where :math:`\omega = \omega_{\boldsymbol{q}j}`, :math:`\omega_{i} = \omega_{\boldsymbol{q}_{i}j_{i}}`, :math:`n = n(\omega)`,
:math:`n_{i} = n(\omega_{i})`, and :math:`V^{(3)}` is the three-phonon matrix element defined in
:ref:`formalism_self_energy`, which vanishes unless :math:`\boldsymbol{q}+\boldsymbol{q}_{1}+\boldsymbol{q}_{2}` is a
reciprocal lattice vector. The numerator of the second term vanishes at :math:`\omega = \omega_{1}+\omega_{2}`, so
the expression has no pole there. Modes with zero frequency (the acoustic modes at :math:`\Gamma`) are excluded from the sum.
In ``MODE = phonons``, the harmonic frequencies and eigenvectors are used at all temperatures. In ``MODE = SCPH``,
the SCPH frequencies :math:`\Omega_{\boldsymbol{q}j}` and eigenvectors of each temperature are used instead, both in
:math:`V^{(3)}` and in the occupation factors, while the cubic force constants are unchanged.
The occupation factors expand as :math:`(1+n)(1+n_{1})(1+n_{2}) - n n_{1} n_{2} = 1 + n + n_{1} + n_{2} + n n_{1} + n_{1} n_{2} + n_{2} n`
and :math:`n(1+n_{1})(1+n_{2}) - (1+n) n_{1} n_{2} = n + n n_{1} + n n_{2} - n_{1} n_{2}`. When ``CLASSICAL = 1``, only the
terms quadratic in the occupations are kept, with :math:`n_{i} = kT/\hbar\omega_{i}`, and the bracket becomes
:math:`4(kT/\hbar)^{2}/(\omega\omega_{1}\omega_{2})`, so that :math:`F^{\mathrm{B}}` is proportional to :math:`T^{2}`
in ``MODE = phonons``.
The bubble contributions to :math:`C_{\mathrm{v}}` and :math:`S` are not computed; they can be obtained by
numerical temperature derivatives of the free energy.


Mean square displacement
------------------------


The displacement-displacement correlation function is given by

.. math::
    :nowrap:

    \begin{align}
     \left< u_{\mu}(0\kappa)u_{\nu}(\ell'\kappa') \right> & = \frac{1}{\sqrt{M_{\kappa}M_{\kappa'}}N_{q}}\sum_{\boldsymbol{q}j}
     \frac{\hbar (2n_{\boldsymbol{q}j}+1)}{2\omega_{\boldsymbol{q}j}}
     \mathrm{Re}\bigg[e_{\mu}(\kappa;\boldsymbol{q}j)e_{\nu}^{*}(\kappa';\boldsymbol{q}j) e^{-i\boldsymbol{q}\cdot\boldsymbol{r}(\ell')}  \bigg],
    \end{align}

where :math:`n_{\boldsymbol{q}j} = 1/(e^{\hbar\omega_{\boldsymbol{q}j}/kT}-1)` is the Bose-Einstein distribution function. When ``UCORR = 1``, the code prints the above correlation function in ``PREFIX``.ucorr. The vector :math:`\boldsymbol{r}(\ell')` is the lattice translation vector to the :math:`\ell'` \ th unit cell, which can be specified by the ``SHIFT_UCORR`` tag.
When ``PRINTMSD`` is turned on, the code prints the mean square displacements, which are the diagonal components of the correlation function:

.. math::

    \left< u_{\mu}^{2}(0\kappa)\right> = \frac{\hbar}{M_{\kappa}N_{q}}\sum_{\boldsymbol{q},j}\frac{1}{\omega_{\boldsymbol{q}j}} |e_{\mu}(\kappa;\boldsymbol{q}j)|^{2}
    \left(n_{\boldsymbol{q}j}+\frac{1}{2}\right).


Phonon DOS
----------

When *KPMODE* = 2, the program *anphon* saves the (one) phonon density of states (DOS) to the file ``PREFIX``.dos.
The one-phonon DOS is given by

.. math::

    \mathrm{DOS}(\omega) = \frac{1}{N_{q}}\sum_{\boldsymbol{q},j}\delta(\omega - \omega_{\boldsymbol{q}j}).

If ``PDOS = 1`` is given, the program also prints the atom-projected phonon DOS which is given by

.. math::
 
    \mathrm{PDOS}(\kappa;\omega) = \frac{1}{N_{q}}\sum_{\boldsymbol{q},j}|\boldsymbol{e}(\kappa;\boldsymbol{q}j)|^{2}\delta(\omega - \omega_{\boldsymbol{q}j}).

In addition, ``TDOS``-tag is available to compute the two-phonon DOS defined by

.. math::

    \mathrm{DOS2}(\omega;\boldsymbol{q};\pm) = \frac{1}{N_{q}}\sum_{\boldsymbol{q}_{1},\boldsymbol{q}_{2}, j_{1}, j_{2}}
    \delta(\omega\pm\omega_{\boldsymbol{q}_{1}j_{1}}-\omega_{\boldsymbol{q}_{2}j_{2}})\delta_{\boldsymbol{q}\pm\boldsymbol{q}_{1},\boldsymbol{q}_{2}+\boldsymbol{G}},

where :math:`\boldsymbol{G}` is a reciprocal lattice vector. The signs :math:`\pm` correspond to absorption and emission processes, respectively. Please note that the computation of the two-phonon DOS can be expensive
especially when :math:`N_{q}` or :math:`N_{\kappa}` is large.


(Atomic) participation ratio
----------------------------

Participation ratio (PR) and atomic participation ratio (APR) defined in the following may be useful to analyze the localized nature of the phonon mode :math:`\boldsymbol{q}j`.

* Participation ratio (PR)

.. math::

    PR_{\boldsymbol{q}j} = \left(\sum_{\kappa}^{N_{\kappa}} \frac{|\boldsymbol{e}(\kappa;\boldsymbol{q}j)|^{2}}{M_{\kappa}}\right)^{2} \Bigg/
    N_{\kappa} \sum_{\kappa}^{N_{\kappa}} \frac{|\boldsymbol{e}(\kappa;\boldsymbol{q}j)|^{4}}{M_{\kappa}^{2}}

* Atomic participation ratio (APR)

.. math::

    APR_{\boldsymbol{q}j,\kappa} = \frac{|\boldsymbol{e}(\kappa;\boldsymbol{q}j)|^{2}}{M_{\kappa}} \Bigg/ \left(  N_{\kappa} \sum_{\kappa}^{N_{\kappa}} \frac{|\boldsymbol{e}(\kappa;\boldsymbol{q}j)|^{4}}{M_{\kappa}^{2}} \right)^{1/2}

For an extended eigenmode, the PR value is of order 1, whereas for a localized eigenmodes PR is of order :math:`1/N_{\kappa}` [Hafner1993]_. APR is an atomic decomposition of PR that satisfies :math:`PR_{\boldsymbol{q}j} = \sum_{\kappa} (APR_{\boldsymbol{q}j,\kappa})^{2}`. To print the PR and APR, please set ``MODE = phonons`` and ``PRINTPR = 1`` in the ``&analysis`` entry field. 

Scattering phase space
-----------------------

When *KPMODE* = 2 and ``SPS = 1``, the three-phonon scattering phase space :math:`P_{3}` is calculated and saved to the file ``PREFIX``.sps. :math:`P_{3}` is defined as

.. math::
    
    P_{3}(\boldsymbol{q}j) = \frac{1}{3m^{3}} (2P_{3}^{(+)}(\boldsymbol{q}j) + P_{3}^{(-)}(\boldsymbol{q}j)),

where :math:`m` is the number of phonon branches and 

.. math::
    
    P_{3}^{(\pm)}(\boldsymbol{q}j) = \frac{1}{N_{q}}\sum_{\boldsymbol{q}_{1},\boldsymbol{q}_{2}, j_{1}, j_{2}}\delta(\omega_{\boldsymbol{q}j}\pm\omega_{\boldsymbol{q}_{1}j_{1}}-\omega_{\boldsymbol{q}_{2}j_{2}})\delta_{\boldsymbol{q}\pm\boldsymbol{q}_{1},\boldsymbol{q}_{2}+\boldsymbol{G}}.

*anphon* also prints the total scattering phase space

.. math::

    P_{3} = \frac{1}{N_{q}}\sum_{\boldsymbol{q}j} P_{3}(\boldsymbol{q}j).

When ``SPS = 2``, the three-phonon scattering phase space with the occupation factor :math:`W_{3}^{(\pm)}` will be calculated and saved to the file ``PREFIX``.sps_Bose. :math:`W_{3}^{(\pm)}` is defined as

.. math::

    W_{3}^{(\pm)}(\boldsymbol{q}j) = \frac{1}{N_{q}}{\sum_{\boldsymbol{q}_{1},\boldsymbol{q}_{2}, j_{1}, j_{2}}}
    \left\{
      \begin{array}{c}
      n_{2} - n_{1} \\
      n_{1} + n_{2} + 1
      \end{array}
    \right\}
    \delta(\omega_{\boldsymbol{q}j}\pm\omega_{\boldsymbol{q}_{1}j_{1}}-\omega_{\boldsymbol{q}_{2}j_{2}})\delta_{\boldsymbol{q}\pm\boldsymbol{q}_{1},\boldsymbol{q}_{2}+\boldsymbol{G}}.

Here, :math:`n_{1}=n(\omega_{\boldsymbol{q}_{1}j_{1}})` and :math:`n_{2}=n(\omega_{\boldsymbol{q}_{2}j_{2}})` where :math:`n(\omega) = \frac{1}{e^{\hbar\omega/k_B T}-1}` is the Bose-Einstein distribution function. Since :math:`n(\omega)` is temperature dependent, :math:`W_{3}^{(\pm)}` is also temperature dependent. The file ``PREFIX``.sps_Bose contains :math:`W_{3}^{(\pm)}` for all phonon modes at various temperatures specified with ``TMIN``, ``TMAX``, and ``DT`` tags.

.. _formalism_gruneisen:

Gr\ |umulaut_u|\ neisen parameter
---------------------------------

The volumetric mode Gr\ |umulaut_u|\ neisen parameter, defined as :math:`\gamma_{\boldsymbol{q}j} = - \frac{\partial \log{\omega_{\boldsymbol{q}j}}}{\partial \log{V}}`,
is calculated by

.. math::

    \gamma_{\boldsymbol{q}j}= -\frac{(\boldsymbol{e}_{\boldsymbol{q}j}^{*})^{\mathrm{T}} \delta D(\boldsymbol{q})\boldsymbol{e}_{\boldsymbol{q}j}}{6\omega_{\boldsymbol{q}j}^{2}},

where :math:`\delta D(\boldsymbol{q})` is a change in the dynamical matrix due to a volume change :math:`\delta V`,
which is given by

.. math::
    :nowrap:

    \begin{align}
     \delta D_{\mu\nu}(\kappa\kappa^{\prime};\boldsymbol{q}) &= \frac{1}{\sqrt{M_{\kappa}M_{\kappa^{\prime}}}}
     \sum_{\ell^{\prime}}\delta\Phi_{\mu\nu}(\ell\kappa;\ell^{\prime}\kappa^{\prime})\exp{\left[i\boldsymbol{q}\cdot(\boldsymbol{r}(\ell^{\prime})-\boldsymbol{r}(\ell))\right]},\\
     \delta\Phi_{\mu\nu}(\ell\kappa;\ell^{\prime}\kappa^{\prime})
     &= \sum_{\ell^{\prime\prime},\kappa^{\prime\prime},\lambda}\Phi_{\mu\nu\lambda}(\ell\kappa;\ell^{\prime}\kappa^{\prime};\ell^{\prime\prime}\kappa^{\prime\prime})r_{\lambda}(\ell^{\prime\prime}\kappa^{\prime\prime}).
    \end{align}

Please set ``GRUNEISEN = 1`` and give an appropriate ``FCSFILE`` file containing cubic IFCs to print the volumetric Gr\ |umulaut_u|\ neisen parameters.

The **generalized** Gr\ |umulaut_u|\ neisen parameters resolve the frequency change with respect to each component of the symmetric strain tensor :math:`\varepsilon_{\mu\nu}`:

.. math::

    \gamma_{\boldsymbol{q}j}^{\mu\nu} = - \frac{\partial \log{\omega_{\boldsymbol{q}j}}}{\partial \varepsilon_{\mu\nu}}
    = -\frac{1}{2\omega_{\boldsymbol{q}j}^{2}} (\boldsymbol{e}_{\boldsymbol{q}j}^{*})^{\mathrm{T}} \frac{\partial D(\boldsymbol{q})}{\partial \varepsilon_{\mu\nu}}\boldsymbol{e}_{\boldsymbol{q}j},

where the strain derivative of the harmonic IFCs is obtained from the cubic IFCs as

.. math::

    \frac{\partial \Phi_{\alpha\beta}(\ell\kappa;\ell^{\prime}\kappa^{\prime})}{\partial \varepsilon_{\mu\nu}}
    = \sum_{\ell^{\prime\prime},\kappa^{\prime\prime}}\Phi_{\alpha\beta\mu}(\ell\kappa;\ell^{\prime}\kappa^{\prime};\ell^{\prime\prime}\kappa^{\prime\prime})r_{\nu}(\ell^{\prime\prime}\kappa^{\prime\prime}),

and the off-diagonal components are symmetrized as :math:`(\gamma_{\boldsymbol{q}j}^{\mu\nu} + \gamma_{\boldsymbol{q}j}^{\nu\mu})/2`.
The volumetric parameter is recovered as one third of the trace,
:math:`\gamma_{\boldsymbol{q}j} = (\gamma_{\boldsymbol{q}j}^{xx} + \gamma_{\boldsymbol{q}j}^{yy} + \gamma_{\boldsymbol{q}j}^{zz})/3`.
The generalized parameters are useful for studying anisotropic thermal expansion and thermal stress within the quasi-harmonic approximation.
Please set ``GRUNEISEN = 2`` (diagonal components) or ``GRUNEISEN = 3`` (all 6 components in the Voigt order xx, yy, zz, yz, xz, xy) to print them.


.. _formalism_self_energy:

Anharmonic self-energy
-----------------------

The anharmonic self-energy due to cubic anharmonicity to the lowest order is given by

.. math::
    :label: self3

    \Sigma_{\boldsymbol{q}j}(i\omega_m) &= \frac{1}{2\hbar^{2}}\sum_{\boldsymbol{q}_{1},\boldsymbol{q}_{2}}\sum_{j_{1},j_{2}}
    |V^{(3)}_{-\boldsymbol{q}j,\boldsymbol{q}_{1}j_{1},\boldsymbol{q}_{2}j_{2}}|^{2} \notag \\
    & \times \left[ \frac{n_{1}+n_{2} + 1}{i\omega_{m} + \omega_{1} + \omega_{2}} - \frac{n_{1}+n_{2} + 1}{i\omega_{m} - \omega_{1} - \omega_{2}} 
    + \frac{n_{1}-n_{2}}{i\omega_{m} - \omega_{1} + \omega_{2}} - \frac{n_{1}-n_{2}}{i\omega_{m} + \omega_{1} - \omega_{2}} \right],
    
where :math:`i\omega_{m}` is the Matsubara frequency. In equation :eq:`self3`, we simply denoted :math:`\omega_{\boldsymbol{q}_{i}j_{i}}` as :math:`\omega_{i}`. The matrix element :math:`V^{(3)}` is given by

.. math::
  
    V^{(3)}_{\boldsymbol{q}j,\boldsymbol{q}^{\prime}j^{\prime},\boldsymbol{q}^{\prime\prime}j^{\prime\prime}} 
    & = \left( \frac{\hbar}{2N_{q}}\right)^{\frac{3}{2}}
    \frac{1}{\sqrt{\omega_{\boldsymbol{q}j}\omega_{\boldsymbol{q}^{\prime}j^{\prime}}\omega_{\boldsymbol{q}^{\prime\prime}j^{\prime\prime}}}}
    \sum_{\ell,\ell^{\prime},\ell^{\prime\prime}}
    \exp{\left[\mathrm{i}(\boldsymbol{q}\cdot\boldsymbol{r}(\ell)+\boldsymbol{q}^{\prime}\cdot\boldsymbol{r}(\ell^{\prime})+\boldsymbol{q}^{\prime\prime}\cdot\boldsymbol{r}(\ell^{\prime\prime}))\right]} \notag \\
    & \times \sum_{\kappa,\kappa^{\prime},\kappa^{\prime\prime}} \frac{1}{\sqrt{M_{\kappa}M_{\kappa^{\prime}}M_{\kappa^{\prime\prime}}}}
    \sum_{\mu,\nu,\lambda}
    \Phi_{\mu\nu\lambda}(\ell\kappa;\ell^{\prime}\kappa^{\prime};\ell^{\prime\prime}\kappa^{\prime\prime}) 
    e_{\mu}(\kappa;\boldsymbol{q}j)e_{\nu}(\kappa^{\prime};\boldsymbol{q}^{\prime}j^{\prime})e_{\lambda}(\kappa^{\prime\prime};\boldsymbol{q}^{\prime\prime}j^{\prime\prime}) \; ,
    

which becomes zero unless :math:`\boldsymbol{q}+\boldsymbol{q}^{\prime}+\boldsymbol{q}^{\prime\prime}` is an integral multiple of :math:`\boldsymbol{G}=n_{1}\boldsymbol{b}_{1}+n_{2}\boldsymbol{b}_{2}+n_{3}\boldsymbol{b}_{3}`.
Phonon linewidth :math:`\Gamma_{\boldsymbol{q}j}`, which is the imaginary part of the phonon self-energy, can be obtained by the analytic continuation to the real axis (:math:`i\omega_{m}\to \omega + i0^{+}`) as

.. math::
    :label: selfmod

     \Gamma^{\mathrm{3rd}}_{\boldsymbol{q}j}(\omega) &= \frac{\pi}{2\hbar^{2}}\sum_{\boldsymbol{q}_{1},\boldsymbol{q}_{2}}\sum_{j_{1},j_{2}}
     |V^{(3)}_{-\boldsymbol{q}j,\boldsymbol{q}_{1}j_{1},\boldsymbol{q}_{2}j_{2}}|^{2} \notag \\
     & \times \left[ -(n_{1}+n_{2} + 1)\delta{(\omega + \omega_{1} + \omega_{2})} + (n_{1}+n_{2} + 1) \delta{(\omega - \omega_{1} - \omega_{2})} \right. \notag \\
     & \left. \hspace{12mm} - (n_{1}-n_{2})\delta{(\omega - \omega_{1} + \omega_{2})} + (n_{1}-n_{2})\delta{(\omega + \omega_{1} - \omega_{2})} \right].

The computation of equation :eq:`selfmod` is the most expensive part of the thermal conductivity calculations.
Therefore, we employ the crystal symmetry to reduce the number of triplet pairs :math:`(\boldsymbol{q}j,\boldsymbol{q}^{\prime}j^{\prime},\boldsymbol{q}^{\prime\prime}j^{\prime\prime})` of :math:`V^{(3)}` to calculate.
The tetrahedron integration weights (``ISMEAR = -1``) are not invariant under the point group; they are therefore averaged
over each group of symmetry-equivalent triplets, and the reduced calculation agrees closely with the unreduced one.
The same reduction is used for the collision operator of the iterative BTE solvers.
To disable the reduction, please set ``TRISYM = 0``.

.. _formalism_4ph:

Four-phonon scattering
~~~~~~~~~~~~~~~~~~~~~~

When ``INCLUDE_4PH = 1`` in the ``&kappa`` field, the scattering rates due to the quartic anharmonicity are added
to those of equation :eq:`selfmod` [Feng2016]_. They come from the imaginary part of the lowest-order self-energy
diagram with two quartic vertices (the "sunset" diagram), whose matrix element is defined in the same way as
:math:`V^{(3)}`:

.. math::

    V^{(4)}_{\boldsymbol{q}j,\boldsymbol{q}_{1}j_{1},\boldsymbol{q}_{2}j_{2},\boldsymbol{q}_{3}j_{3}}
    & = \left( \frac{\hbar}{2N_{q}}\right)^{2}
    \frac{1}{\sqrt{\omega_{\boldsymbol{q}j}\omega_{1}\omega_{2}\omega_{3}}}
    \sum_{\ell,\ell_{1},\ell_{2},\ell_{3}}
    \exp{\left[\mathrm{i}(\boldsymbol{q}\cdot\boldsymbol{r}(\ell)+\boldsymbol{q}_{1}\cdot\boldsymbol{r}(\ell_{1})+\boldsymbol{q}_{2}\cdot\boldsymbol{r}(\ell_{2})+\boldsymbol{q}_{3}\cdot\boldsymbol{r}(\ell_{3}))\right]} \notag \\
    & \times \sum_{\kappa,\kappa_{1},\kappa_{2},\kappa_{3}} \frac{1}{\sqrt{M_{\kappa}M_{\kappa_{1}}M_{\kappa_{2}}M_{\kappa_{3}}}}
    \sum_{\mu,\nu,\lambda,\rho}
    \Phi_{\mu\nu\lambda\rho}(\ell\kappa;\ell_{1}\kappa_{1};\ell_{2}\kappa_{2};\ell_{3}\kappa_{3})
    e_{\mu}(\kappa;\boldsymbol{q}j)e_{\nu}(\kappa_{1};\boldsymbol{q}_{1}j_{1})e_{\lambda}(\kappa_{2};\boldsymbol{q}_{2}j_{2})e_{\rho}(\kappa_{3};\boldsymbol{q}_{3}j_{3}),

which vanishes unless :math:`\boldsymbol{q}+\boldsymbol{q}_{1}+\boldsymbol{q}_{2}+\boldsymbol{q}_{3}` is a reciprocal lattice vector.
The four-phonon linewidth is then

.. math::
    :label: gamma4ph

     \Gamma^{\mathrm{4ph}}_{\boldsymbol{q}j}(\omega) &= \frac{\pi}{6\hbar^{2}}\sum_{\boldsymbol{q}_{1},\boldsymbol{q}_{2},\boldsymbol{q}_{3}}\sum_{j_{1},j_{2},j_{3}}
     |V^{(4)}_{-\boldsymbol{q}j,\boldsymbol{q}_{1}j_{1},\boldsymbol{q}_{2}j_{2},\boldsymbol{q}_{3}j_{3}}|^{2} \notag \\
     & \times \Big\{ \left[(1+n_{1})(1+n_{2})(1+n_{3}) - n_{1}n_{2}n_{3}\right]\delta(\omega-\omega_{1}-\omega_{2}-\omega_{3}) \notag \\
     & \hspace{5mm} + \left[n_{1}(1+n_{2})(1+n_{3}) - (1+n_{1})n_{2}n_{3}\right]\left[\delta(\omega+\omega_{1}-\omega_{2}-\omega_{3}) - \delta(\omega-\omega_{1}+\omega_{2}+\omega_{3})\right] \notag \\
     & \hspace{5mm} + \left[n_{2}(1+n_{3})(1+n_{1}) - (1+n_{2})n_{3}n_{1}\right]\left[\delta(\omega+\omega_{2}-\omega_{3}-\omega_{1}) - \delta(\omega-\omega_{2}+\omega_{3}+\omega_{1})\right] \notag \\
     & \hspace{5mm} + \left[n_{3}(1+n_{1})(1+n_{2}) - (1+n_{3})n_{1}n_{2}\right]\left[\delta(\omega+\omega_{3}-\omega_{1}-\omega_{2}) - \delta(\omega-\omega_{3}+\omega_{1}+\omega_{2})\right] \Big\},

where :math:`n_{i} = n(\omega_{i})`. The first line describes the splitting of the phonon :math:`\boldsymbol{q}j` into three phonons
(:math:`1\to 3`). In each of the remaining three lines, the first delta function describes the processes in which the phonon
:math:`\boldsymbol{q}j` and one phonon of the quartet are annihilated and the other two are created (:math:`2\leftrightarrow 2`),
and the second one the processes in which the phonon :math:`\boldsymbol{q}j` and two phonons of the quartet merge into the
third one (:math:`3\to 1`).
The occupation factors expand as :math:`(1+n_{1})(1+n_{2})(1+n_{3}) - n_{1}n_{2}n_{3} = 1+n_{1}+n_{2}+n_{3}+n_{1}n_{2}+n_{2}n_{3}+n_{3}n_{1}`
and :math:`n_{1}(1+n_{2})(1+n_{3}) - (1+n_{1})n_{2}n_{3} = n_{1}+n_{1}n_{2}+n_{1}n_{3}-n_{2}n_{3}`; the leading 1 is dropped when ``CLASSICAL = 1``.
The partner term :math:`-\delta(\omega+\omega_{1}+\omega_{2}+\omega_{3})` of the first line, which makes
:math:`\Gamma^{\mathrm{4ph}}(\omega)` odd in :math:`\omega`, never contributes for :math:`\omega > 0` and is omitted.
In the thermal conductivity calculation, :math:`\Gamma^{\mathrm{4ph}}_{\boldsymbol{q}j} = \Gamma^{\mathrm{4ph}}_{\boldsymbol{q}j}(\omega_{\boldsymbol{q}j})`
is added to the out-scattering rate only (see :ref:`below <kappa_peierls>`).

Because the number of quartets grows as :math:`N_{q}^{2}` per mode instead of :math:`N_{q}` for the triplets,
equation :eq:`gamma4ph` is usually evaluated on a :math:`\boldsymbol{q}` mesh coarser than that of the three-phonon
part (``KMESH_COARSE``), and the resulting linewidths are interpolated onto the dense mesh by the method selected
with ``INTERPOLATOR``. The delta functions are evaluated with the Gaussian, Lorentzian, or adaptive smearing chosen
with ``ISMEAR_4PH`` and ``EPSILON_4PH`` (see :ref:`formalism_delta_function`); the tetrahedron method is not available.
For each quartet the band triples :math:`(j_{1},j_{2},j_{3})` outside the smearing window are discarded first, and
the matrix elements of the remaining ones are obtained by successive contractions of the Fourier-transformed
quartic force constants with the eigenvectors (see :ref:`INCLUDE_4PH <anphon_include_4ph>`). Quartets equivalent under the
little group of :math:`\boldsymbol{q}` or under the permutations of :math:`(\boldsymbol{q}_{1},\boldsymbol{q}_{2},\boldsymbol{q}_{3})`
are evaluated only once. The four-phonon linewidths are written to ``PREFIX.kappa.h5`` (``/scattering/4ph/gamma``).


Isotope scattering
------------------

The effect of isotope scatterings can be considered by the mass perturbation approach proposed by S. Tamura [Tamura1983]_ by the ``ISOTOPE``-tag.
The corresponding phonon linewidth is given by

.. math::

    \Gamma_{\boldsymbol{q}j}^{\mathrm{iso}}(\omega)= \frac{\pi}{4N_{q}} \omega_{\boldsymbol{q}j}^{2}\sum_{\boldsymbol{q}_{1},j_{1}}\delta(\omega-\omega_{\boldsymbol{q}_{1}j_{1}})
    \sum_{\kappa}g_{2}(\kappa)|\boldsymbol{e}^{*}(\kappa;\boldsymbol{q}_{1}j_{1})\cdot\boldsymbol{e}(\kappa;\boldsymbol{q}j)|^{2},

where :math:`g_{2}` is a dimensionless factor given by

.. math::

    g_{2}(\kappa)=\sum_{i}f_{i}(\kappa)\left(1 - \frac{m_{i}(\kappa)}{M_{\kappa}}\right)^{2}.

Here, :math:`f_{i}` is the fraction of the :math:`i`\ th isotope of an element having mass :math:`m_i`, 
and :math:`M_{\kappa}=\sum_{i}f_{i}m_{i}(\kappa)` is the average mass, respectively.
The :math:`g_{2}` values should be provided by the ``ISOFACT``-tag.
The average mass :math:`M_{\kappa}` is substituted by the value specified in the ``MASS``-tag.

.. _kappa:

Lattice thermal conductivity
----------------------------

The lattice thermal conductivity computed by *anphon* is the sum of a particle-like (Peierls) term, in which heat is carried
by the populations of the individual phonon modes, and, when requested, a wave-like (coherent) term arising from the
interband elements of the heat-flux operator [Simoncelli2022]_. The Peierls term can be evaluated at several levels of approximation to
the phonon Boltzmann transport equation, selected by the ``SOLVER`` tag of the ``&kappa`` field.

.. _kappa_peierls:

Peierls term
~~~~~~~~~~~~

Under a small temperature gradient the steady-state phonon distribution is written as

.. math::

    n_{\lambda} = n^{0}_{\lambda} - \frac{\partial n^{0}_{\lambda}}{\partial T}\,\boldsymbol{F}_{\lambda}\cdot\nabla T,

where :math:`\lambda=\boldsymbol{q}j`, :math:`n^{0}` is the Bose--Einstein distribution and :math:`\boldsymbol{F}_{\lambda}` is the
*deviation function*, a vector with the dimension of length (the mean free displacement of mode :math:`\lambda`). Linearizing
the collision term of the Boltzmann transport equation (BTE) [Omini1996]_ [Fugallo2013]_ [Chaput2013]_ turns it into a linear system for :math:`\boldsymbol{F}`:

.. math::
    :label: lbte

    \sum_{\lambda'}\Omega_{\lambda\lambda'}\,\boldsymbol{F}_{\lambda'} = \boldsymbol{v}_{\lambda}.

The collision operator :math:`\Omega` has the dimension of a scattering rate. Its diagonal is the out-scattering rate
:math:`\Omega_{\lambda\lambda}=\tau_{\lambda}^{-1}`, and its off-diagonal elements :math:`\Omega_{\lambda\lambda'}`
(:math:`\lambda\neq\lambda'`) are the in-scattering terms, built from the same three-phonon matrix elements
:math:`|V^{(3)}|^{2}` and energy-conserving :math:`\delta` functions as :math:`\Gamma^{\mathrm{3rd}}` in equation :eq:`selfmod`.
Once :math:`\boldsymbol{F}` is known, the Peierls term of the thermal conductivity is

.. math::
    :label: kappa_lbte

    \kappa_{\mathrm{ph}}^{\mu\nu}(T) = \frac{1}{VN_{q}}\sum_{\lambda} c_{\lambda}\,v^{\mu}_{\lambda}F^{\nu}_{\lambda},

where :math:`V` is the unit cell volume and :math:`c_{\lambda} = \hbar\omega_{\lambda}\partial n^{0}_{\lambda}/\partial T`.
Isotope scattering enters :math:`\Omega` with its own in-scattering term when ``ISOTOPE_INSCATTERING = 1`` (default);
boundary scattering (``LEN_BOUNDARY``) and :ref:`four-phonon scattering <formalism_4ph>` (``INCLUDE_4PH = 1``) are added to the diagonal only.
The four solvers below differ in how equation :eq:`lbte` is solved.

Relaxation-time approximation (``SOLVER = RTA``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Dropping the in-scattering part of :math:`\Omega` gives :math:`\boldsymbol{F}_{\lambda}=\boldsymbol{v}_{\lambda}\tau_{\lambda}`, and equation :eq:`kappa_lbte` becomes

.. math::
  
  \kappa_{\mathrm{ph}}^{\mu\nu}(T) = \frac{1}{V N_{q}} \sum_{\boldsymbol{q},j}c_{\boldsymbol{q}j}(T)\,W_{\boldsymbol{q}j}^{\mu\nu}\,\tau_{\boldsymbol{q}j}(T),
  \qquad
  W_{\boldsymbol{q}j}^{\mu\nu} = \sum_{j'\in\mathcal{B}(j)} v_{\boldsymbol{q}jj'}^{\mu}v_{\boldsymbol{q}j'j}^{\nu},

where :math:`V` is the unit cell volume, :math:`c_{\boldsymbol{q}j} = \hbar\omega_{\boldsymbol{q}j}\partial n_{\boldsymbol{q}j}/\partial T`, :math:`\tau_{\boldsymbol{q}j}(T)` is the phonon lifetime,
and :math:`\mathcal{B}(j)` is the set of branches degenerate with :math:`j` at :math:`\boldsymbol{q}` (see below).
For a non-degenerate branch :math:`W_{\boldsymbol{q}j}^{\mu\nu} = v_{\boldsymbol{q}j}^{\mu}v_{\boldsymbol{q}j}^{\nu}` and the familiar expression is recovered.
For a degenerate multiplet the sum of :math:`W` over its members is

.. math::

    \sum_{j\in\mathcal{B}} W_{\boldsymbol{q}j}^{\mu\nu}
    = \sum_{j,j'\in\mathcal{B}} v_{\boldsymbol{q}jj'}^{\mu}v_{\boldsymbol{q}j'j}^{\nu}
    = \mathrm{Tr}\left(P_{\mathcal{B}}V^{\mu}P_{\mathcal{B}}V^{\nu}\right),

where :math:`V^{\mu}` is the matrix with elements :math:`v_{\boldsymbol{q}jj'}^{\mu}` over all branches and
:math:`P_{\mathcal{B}}` is the projector onto the multiplet, i.e. in the band basis the diagonal matrix equal to 1 for
:math:`j\in\mathcal{B}` and 0 otherwise. A unitary rotation :math:`\mathcal{W}` of the eigenvectors inside the multiplet maps
:math:`V^{\mu}\to\mathcal{W}^{\dagger}V^{\mu}\mathcal{W}` and leaves :math:`P_{\mathcal{B}}` unchanged, so the trace is
invariant, whereas the individual products :math:`v_{\boldsymbol{q}jj}^{\mu}v_{\boldsymbol{q}jj}^{\nu}` are not. The pairs :math:`j\neq j'` inside one multiplet are counted here and are
excluded from the :ref:`coherent term <kappa_coherent>`, so nothing is double counted; the two together form the invariant.

This assignment is a convention, not a physical distinction. At exact degeneracy the coherent expression of Ref. [Simoncelli2022]_ reduces
to the Peierls form for such a pair, so the total conductivity is the same whichever term the pair is assigned to; the assignment
only makes each part separately basis invariant. The *particle-like* part defined here therefore equals the *populations* of
Ref. [Simoncelli2022]_ plus the coherences between degenerate branches, and differs from their populations wherever multiplets occur. For
symmetry-enforced degeneracies the multiplet is a well-defined object and the definition is exact; for accidental near-degeneracies
it depends on the numerical tolerance below, whose effect on the total is the factor reported by the warning described there.
Only the total should be compared between calculations or with Ref. [Simoncelli2022]_.

The multiplets are detected from the frequencies only. The frequency-degeneracy groups over which the self-energy is averaged
(tolerance :math:`10^{-7}` Ry) are subdivided with an anchored tolerance of :math:`10^{-6}\ \mathrm{cm}^{-1}`, so that the lifetime is
constant within every multiplet by construction. This is an empirical numerical criterion: a pair whose true splitting
:math:`\Delta\omega` falls below it is treated as degenerate, which overestimates its weight by the factor
:math:`1+[\Delta\omega/(\Gamma_{\boldsymbol{q}j}+\Gamma_{\boldsymbol{q}j'})]^{2}` relative to the coherent expression. *anphon* reports the
largest such ratio among the merged multiplets and warns when it exceeds 0.1.

Because :math:`W` is contracted for each :math:`\boldsymbol{q}` as soon as the velocity matrix at that point is formed, the full matrix
(of size :math:`N_{q}(3N_{\kappa})^{2}\times3`) is never stored unless the coherent term is requested.

The per-mode linewidths written to ``PREFIX.kappa.h5`` (``/scattering/3ph/gamma``, ``/scattering/4ph/gamma``) and to the
text ``.result`` files are averaged over the same frequency-degeneracy groups before they are written, as the isotope linewidths
are. A single mode of a degenerate multiplet has no basis-independent linewidth: rotating the eigenvectors inside the
multiplet (a perturbation of :math:`10^{-7}` of the force constants suffices) changes individual values by up to 80 %, while their
mean is unchanged. Since the RTA conductivity uses the averaged values anyway, it is unaffected; the four-phonon rates that the
iterative solvers add to the diagonal of the collision operator are averaged in the same way. The written values are complete
for each irreducible :math:`\boldsymbol{q}` only, and a restart keeps only complete :math:`\boldsymbol{q}`: the modes of a
partially written one (in an interrupted file, or a text ``.result`` file of version 1.x) are recomputed. The ``gamma`` datasets carry the
attribute ``averaging = degenerate_modes_mean_1e-7Ry``; raw per-mode rates are not stored, so the ``--noavg`` option of
``tools/analyzer.py`` has an effect only on ``.result`` files of version 1.x.

The phonon lifetime is estimated using the Matthiessen's rule as

.. math::

    \tau_{\boldsymbol{q}j}^{-1}(T) = 2 (\Gamma_{\boldsymbol{q}j}^{\mathrm{3rd}}(T) + \Gamma_{\boldsymbol{q}j}^{\mathrm{4ph}}(T) + \Gamma_{\boldsymbol{q}j}^{\mathrm{iso}}),

where the four-phonon term :eq:`gamma4ph` is included only when ``INCLUDE_4PH = 1``.

When boundary scattering is included (``LEN_BOUNDARY``), the additional rate :math:`|\boldsymbol{v}_{\boldsymbol{q}j}|/L` uses, for a
degenerate multiplet, the block mean speed :math:`|\boldsymbol{v}|^{2}=d_{\mathcal{B}}^{-1}\sum_{j\in\mathcal{B}}\sum_{\mu}W_{\boldsymbol{q}j}^{\mu\mu}`,
which is likewise basis invariant.

.. note::

   With this formulation the conductivity does not depend on the choice of unit cell: a primitive cell and a commensurate
   supercell sampling the same :math:`\boldsymbol{q}` points give identical :math:`\kappa` (to :math:`10^{-15}` relative in a
   silicon test) when a fixed-width smearing (``ISMEAR = 0`` or ``1``) is used. The tetrahedron method (``ISMEAR = -1``)
   interpolates frequencies at fixed sorted branch indices and its linewidths remain cell dependent at the few-percent level; the
   adaptive smearing (``ISMEAR = 2``) and the iterative solvers (``SOLVER = IBTE`` and variants) still use finite-difference
   velocities and are cell dependent as well.

The lattice thermal conductivity is saved in the file ``PREFIX``.kl.

The spectra of the lattice thermal conductivity :math:`\kappa_{\mathrm{ph}}^{\mu\mu}(\omega)` can also be calculated by setting ``KAPPA_SPEC = 1`` in the ``&kappa`` field. :math:`\kappa_{\mathrm{ph}}^{\mu\mu}(\omega)` is defined as 

.. math::
    \kappa_{\mathrm{ph}}^{\mu\mu}(\omega) = \frac{1}{V N_{q}}\sum_{\boldsymbol{q},j}c_{\boldsymbol{q}j}v_{\boldsymbol{q}j}^{\mu}v_{\boldsymbol{q}j}^{\mu}\tau_{\boldsymbol{q}j} \delta(\omega-\omega_{\boldsymbol{q}j}).

If we integrate this quantity over :math:`\omega`, we then obtain the bulk thermal conductivity, namely :math:`\kappa_{\mathrm{ph}}^{\mu\mu} = \int_{0}^{\infty} \kappa_{\mathrm{ph}}^{\mu\mu}(\omega) \; \mathrm{d}\omega`.

The RTA treats every scattering event as if it simply destroyed the excess of phonons in mode :math:`\lambda`. In reality a
collision also *creates* excess phonons in the partner modes, and for momentum-conserving (Normal) processes this in-scattering
keeps the heat current flowing. Neglecting it makes the RTA underestimate :math:`\kappa`, sometimes by a large factor in
materials with weak Umklapp scattering (diamond, BAs, graphene, or any crystal at low temperature). The three solvers below keep
the in-scattering term of :math:`\Omega`; for the same mesh and smearing they converge to the same :math:`\kappa`, which is
written to ``PREFIX``.kl_iter.

.. _kappa_beyond_rta:

Iterative solution (``SOLVER = IBTE``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Starting from the RTA, :math:`\boldsymbol{F}^{(0)}_{\lambda}=\boldsymbol{v}_{\lambda}\tau_{\lambda}`,
the in-scattering term is evaluated with the previous iterate and the diagonal is inverted:

.. math::

    \boldsymbol{F}^{(i+1)}_{\lambda} = \tau_{\lambda}\Big[\boldsymbol{v}_{\lambda} - \sum_{\lambda'\neq\lambda}\Omega_{\lambda\lambda'}\boldsymbol{F}^{(i)}_{\lambda'}\Big].

The new iterate is mixed with the old one, :math:`\boldsymbol{F}\leftarrow\alpha\boldsymbol{F}^{(i+1)}+(1-\alpha)\boldsymbol{F}^{(i)}`
with :math:`\alpha` = ``IBTE_MIXING``, and the iteration stops when the relative residual
:math:`\|\Omega\boldsymbol{F}-\boldsymbol{v}\|/\|\boldsymbol{v}\|` falls below ``ITER_THRESHOLD`` (after at least
``MIN_CYCLE`` and at most ``MAX_CYCLE`` iterations). The residual of the current iterate comes at no extra cost from the
update itself; the stored :math:`\boldsymbol{F}` is the iterate whose residual was tested, not the update computed
alongside it.
This scheme is simple and each iteration costs one application of :math:`\Omega`, but it converges slowly, or oscillates,
when in-scattering is strong (low temperature, very high :math:`\kappa`). A mixing factor below 1 damps the oscillation.

Variational solution (``SOLVER = VBTE``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because :math:`\Omega` is symmetric and positive semidefinite (a consequence of detailed
balance), the solution of equation :eq:`lbte` is also the extremum of the quadratic functional
:math:`\mathcal{F}[\boldsymbol{F}]=\tfrac{1}{2}\boldsymbol{F}\cdot\Omega\boldsymbol{F}-\boldsymbol{F}\cdot\boldsymbol{v}`, and :math:`\kappa` is
its value at the extremum [Fugallo2013]_. *anphon* finds it with the preconditioned conjugate-gradient method. Two practical
consequences follow: the iteration cannot oscillate, and an approximate :math:`\boldsymbol{F}` with residual
:math:`\epsilon=\|\Omega\boldsymbol{F}-\boldsymbol{v}\|/\|\boldsymbol{v}\|` gives :math:`\kappa` with an error of order :math:`\epsilon^{2}`.
For ``VBTE``, ``ITER_THRESHOLD`` is therefore the tolerance on this relative residual, and ``MIN_CYCLE`` and ``IBTE_MIXING``
are ignored. ``VBTE`` usually needs far fewer iterations than ``IBTE`` and is the recommended choice when ``IBTE`` converges
slowly. The tetrahedron and adaptive integration weights make the discretized :math:`\Omega` slightly non-symmetric, and
the conjugate-gradient iteration then stagnates or breaks down at low temperature and when Normal processes dominate.
For these weights (``ISMEAR = -1, 2``), ``VBTE`` therefore uses the BiCGSTAB method, which does not require a symmetric
operator; the variational property is lost, so the error of :math:`\kappa` is then linear in the residual.

Direct solution (``SOLVER = DBTE``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The collision operator restricted to the irreducible wedge is assembled as an explicit dense matrix
and fully diagonalized, :math:`\Omega=\sum_{n}\lambda_{n}\boldsymbol{u}_{n}\boldsymbol{u}_{n}^{\mathsf{T}}`, so that
equation :eq:`lbte` is solved exactly as :math:`\boldsymbol{F}=\sum_{n}\lambda_{n}^{-1}(\boldsymbol{u}_{n}\cdot\boldsymbol{v})\boldsymbol{u}_{n}`
[Chaput2013]_. The eigenvalues :math:`\lambda_{n}` are the scattering rates of collective relaxation modes (printed in
cm\ :sup:`-1`). Energy conservation makes one eigenvalue exactly zero, and near-zero eigenvalues with a large overlap onto the
momentum-drift directions signal that Normal processes dominate, the regime where the RTA fails most. Negative eigenvalues
or a strongly asymmetric matrix indicate that the :math:`\boldsymbol{q}` mesh or the smearing width is not adequate. ``DBTE`` is thus
a diagnostic tool: its memory grows with the square of the number of irreducible modes, so it is intended for small meshes.
The matrix is symmetrized before the diagonalization, so ``DBTE`` is an exact reference for the two iterative solvers only
for fixed-width smearing; with the tetrahedron or adaptive weights it solves the symmetrized kernel, which differs from
the one ``IBTE`` and ``VBTE`` solve by terms that vanish as the mesh is refined.

.. note::

   The three non-RTA solvers currently use finite-difference group velocities and are therefore cell dependent at the level
   discussed in the :ref:`velocity section <group_velocity>`. Please check the convergence of :math:`\kappa` with respect
   to the :math:`\boldsymbol{q}` mesh and the smearing width carefully.

.. _kappa_coherent:

Coherent term
~~~~~~~~~~~~~

The coherent components of lattice thermal conductivity (see Ref. [Simoncelli2019]_), which are associated with the band off-diagonal components of the harmonic heat-flux operator, are calculated as 

.. math::
  
  \kappa_{\mathrm{c}}^{\mu\nu}(T) = \frac{1}{V N_{q}} \sum_{\substack{\boldsymbol{q},jj'\\ j\neq j'}}\frac{c_{\boldsymbol{q}j}\omega_{\boldsymbol{q}j'} + c_{\boldsymbol{q}j'}\omega_{\boldsymbol{q}j}}{\omega_{\boldsymbol{q}j}+ \omega_{\boldsymbol{q}j'}}  v_{\boldsymbol{q}jj'}^{\mu}v_{\boldsymbol{q}j'j}^{\nu} \frac{\Gamma_{\boldsymbol{q}j}+\Gamma_{\boldsymbol{q}j'}}{(\omega_{\boldsymbol{q}j}-\omega_{\boldsymbol{q}j'})^{2}+(\Gamma_{\boldsymbol{q}j}+\Gamma_{\boldsymbol{q}j'})^2},

where :math:`c_{\boldsymbol{q}j} = \hbar\omega_{\boldsymbol{q}j}\partial n_{\boldsymbol{q}j}/\partial T` and :math:`\Gamma_{\boldsymbol{q}j}` is the total phonon linewidth (half width) of phonon mode :math:`\boldsymbol{q}j`. 

:math:`\boldsymbol{v}_{\boldsymbol{q}jj'}` is the band off-diagonal velocity matrix :eq:`velmat`, built from
:math:`\tilde{D}` [Allen1993]_ [Simoncelli2022]_. The sum runs over pairs belonging to *different* degenerate multiplets;
pairs inside one multiplet are already contained in the :ref:`Peierls term <kappa_peierls>`. For such a pair
:math:`\omega_{\boldsymbol{q}j}=\omega_{\boldsymbol{q}j'}`, the prefactor reduces to :math:`c_{\boldsymbol{q}j}` and the Lorentzian
factor to :math:`1/(\Gamma_{\boldsymbol{q}j}+\Gamma_{\boldsymbol{q}j'})`; since the lifetime is constant within a multiplet
this equals :math:`1/(2\Gamma_{\boldsymbol{q}j})=\tau_{\boldsymbol{q}j}`, i.e. exactly the Peierls weight
:math:`c_{\boldsymbol{q}j}v_{\boldsymbol{q}jj'}^{\mu}v_{\boldsymbol{q}j'j}^{\nu}\tau_{\boldsymbol{q}j}`, so the two terms together
form the basis-invariant block trace. The particle-like/wave-like *split* is therefore basis dependent while their sum is not; only the total
should be compared between calculations. When ``KAPPA_COHERENT = 1 | 2`` the coherent component is calculated and saved in ``PREFIX``.kl_coherent.
When ``KAPPA_COHERENT = 2``, all components of the coherent term before summation are saved in ``PREFIX``.kc_elem. Requesting the coherent
term stores the full velocity matrix, :math:`N_{q}(3N_{\kappa})^{2}\times3` complex numbers on the root process.

.. _cumulative_kappa:

Cumulative thermal conductivity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The accumulative lattice thermal conductivity :math:`\kappa_{\mathrm{ph,acc}}^{\mu\nu}(L)` is defined as

.. math::
  
  \kappa_{\mathrm{ph,acc}}^{\mu\nu}(L) = \frac{1}{V N_{q}} \sum_{\boldsymbol{q},j}c_{\boldsymbol{q}j}v_{\boldsymbol{q}j}^{\mu}v_{\boldsymbol{q}j}^{\nu}\tau_{\boldsymbol{q}j}\Theta (L-|\boldsymbol{v}_{\boldsymbol{q}j}|\tau_{\boldsymbol{q}j}),

where :math:`\Theta(x)` is the step function. This quantity can be calculated by using the script ``analyzer.py`` with ``--calc cumulative`` flag.
With the default velocity-matrix formulation ``analyzer.py`` replaces :math:`v^{\mu}_{\boldsymbol{q}j}v^{\nu}_{\boldsymbol{q}j}`
by the degenerate-block diad :math:`W^{\mu\nu}_{\boldsymbol{q}j}` stored in ``PREFIX.kappa.h5`` (and
:math:`|\boldsymbol{v}_{\boldsymbol{q}j}|` by :math:`\sqrt{\mathrm{tr}\,W_{\boldsymbol{q}j}}`), so that the
large-:math:`L` limit coincides with :math:`\kappa_{\mathrm{ph}}` of the ``anphon`` run; text ``.result`` files,
which carry no velocity matrix, fall back to the finite-difference velocities.
One can also use another definition for the accumulative thermal conductivity:

.. math::
  
  \kappa_{\mathrm{ph,acc}}^{\mu\nu}(L) = \frac{1}{V N_{q}} \sum_{\boldsymbol{q},j}c_{\boldsymbol{q}j}v_{\boldsymbol{q}j}^{\mu}v_{\boldsymbol{q}j}^{\nu}\tau_{\boldsymbol{q}j}\Theta (L-|v_{\boldsymbol{q}j}^{\mu}|\tau_{\boldsymbol{q}j}).

In this case, the contribution to the total thermal conductivity is limited only from phonon modes whose mean-free-path along the :math:`\mu`\ -direction is smaller than :math:`L`.
To calculate this, please use the ``--calc cumulative2`` flag and specify the direction :math:`\mu` by the ``--direction`` option.

.. _formalism_delta_function:

Delta function
--------------

To compute the phonon DOSs and the imaginary part of phonon self-energies,
it is necessary to evaluate the Brillouin-zone integration containing Dirac's delta function.
For that purpose, we provide four options through the ``ISMEAR``-tag. The four-phonon scattering rates use the
separate ``ISMEAR_4PH`` and ``EPSILON_4PH`` tags with the same meaning.

Lorentzian smearing (``ISMEAR = 0``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The delta function is replaced by the Lorentzian function as

.. math::

    \delta(\omega) \approx \frac{1}{\pi}\frac{\epsilon}{\omega^{2}+\epsilon^{2}}.

Gaussian smearing (``ISMEAR = 1``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The delta function is replaced by the Gaussian function as

.. math::

    \delta(\omega) \approx \frac{1}{\sqrt{\pi}\epsilon}\exp{(-\omega^{2}/\epsilon^{2})},

which decays faster than the Lorentzian function.
For both smearing methods, :math:`\epsilon` should be given by the ``EPSILON``-tag, which must be chosen carefully
to avoid any unscientific results. :math:`\epsilon` should be small enough to capture detailed phonon structures
such as phonon DOS or energy conservation surface related to three-phonon processes, but it should be large
enough to avoid unscientific oscillations. Choosing an appropriate value for :math:`\epsilon` is not a trivial task
since it may depend on the phonon structure and the density of :math:`\boldsymbol{q}` points.

Adaptive smearing (``ISMEAR = 2``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Gaussian form above is used, but the width is chosen automatically for every term of the sum from the group
velocities and the mesh spacing, following the adaptive-broadening scheme of Yates *et al.* [Yates2007]_. For a three-phonon process
in which the energy-conservation condition involves :math:`\omega_{\boldsymbol{q}_1 j_1}\pm\omega_{\boldsymbol{q}_2 j_2}`,
the width is

.. math::

    \sigma = \alpha \sqrt{\frac{1}{12}\sum_{u=1}^{3}\left[\left(\boldsymbol{v}_{\boldsymbol{q}_1 j_1}\mp\boldsymbol{v}_{\boldsymbol{q}_2 j_2}\right)\cdot\frac{\boldsymbol{b}_{u}}{N_{u}}\right]^{2}},

where :math:`\boldsymbol{b}_{u}/N_{u}` are the steps of the :math:`\boldsymbol{q}` mesh along the three reciprocal lattice
vectors, i.e. the width is the spread of the argument of the delta function over one mesh cell. For the phonon DOS
only :math:`\boldsymbol{v}_{\boldsymbol{q}j}` enters.

For four-phonon processes (``ISMEAR_4PH = 2``) the same argument is applied to equation :eq:`gamma4ph`. Since
:math:`\boldsymbol{q}_{3}` is fixed by momentum conservation, :math:`\boldsymbol{q}_{3} = \boldsymbol{q}-\boldsymbol{q}_{1}-\boldsymbol{q}_{2}+\boldsymbol{G}`,
the argument :math:`\omega+s_{1}\omega_{\boldsymbol{q}_1 j_1}+s_{2}\omega_{\boldsymbol{q}_2 j_2}+s_{3}\omega_{\boldsymbol{q}_3 j_3}`
(:math:`s_{i}=\pm 1`) of each delta function varies independently with :math:`\boldsymbol{q}_{1}` and :math:`\boldsymbol{q}_{2}`
over a mesh cell, and the two spreads are added in quadrature:

.. math::

    \sigma = \alpha \sqrt{\frac{1}{12}\sum_{u=1}^{3}\left\{\left[\left(s_{1}\boldsymbol{v}_{\boldsymbol{q}_1 j_1}-s_{3}\boldsymbol{v}_{\boldsymbol{q}_3 j_3}\right)\cdot\frac{\boldsymbol{b}_{u}}{N_{u}}\right]^{2}
    + \left[\left(s_{2}\boldsymbol{v}_{\boldsymbol{q}_2 j_2}-s_{3}\boldsymbol{v}_{\boldsymbol{q}_3 j_3}\right)\cdot\frac{\boldsymbol{b}_{u}}{N_{u}}\right]^{2}\right\}}.

For the splitting process :math:`\delta(\omega-\omega_{1}-\omega_{2}-\omega_{3})`, for instance, the velocity differences
are :math:`\boldsymbol{v}_{1}-\boldsymbol{v}_{3}` and :math:`\boldsymbol{v}_{2}-\boldsymbol{v}_{3}`, whereas for
:math:`\delta(\omega+\omega_{3}-\omega_{1}-\omega_{2})` they are :math:`\boldsymbol{v}_{1}+\boldsymbol{v}_{3}` and
:math:`\boldsymbol{v}_{2}+\boldsymbol{v}_{3}`. Each of the four channels of equation :eq:`gamma4ph` therefore has its own width
for a given band triple. Here :math:`N_{u}` and the group velocities refer to the four-phonon mesh (``KMESH_COARSE``),
so the widths are larger on a coarser mesh, as they should be. For both ``ISMEAR_4PH = 1`` and ``2`` the Gaussian is
truncated at :math:`|x| \geq 2\sigma`, where it has dropped to :math:`e^{-4}\approx 2\%` of its peak value; this
cutoff is what allows the band triples outside the window to be discarded before their matrix elements are formed.

For both the three- and four-phonon processes, the prefactor :math:`\alpha` is set by the ``ADAPTIVE_FACTOR``-tag
(default 1), and a lower bound of :math:`2\times10^{-5}` Ry (about 2.2 cm\ :sup:`-1`) is imposed on :math:`\sigma`.
``EPSILON`` is not used, so no manual convergence test of the width is needed; the result converges with the
:math:`\boldsymbol{q}` mesh alone. The widths are built from finite-difference velocities, so the result retains a weak
dependence on the choice of unit cell (see the :ref:`velocity section <group_velocity>`).

Tetrahedron method (``ISMEAR = -1``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The program *anphon* employs the tetrahedron method [Blochl1994]_ by default
for numerical evaluations of Brillouin zone integration containing :math:`\delta(\omega)`.
When the tetrahedron method is used, the ``EPSILON``-tag is neglected.
We recommend using the tetrahedron method whenever possible. It is not available for the four-phonon scattering rates,
for which ``ISMEAR_4PH = -1`` is automatically switched to the adaptive smearing.

.. _formalism_SCPH:

Self-consistent phonon (SCPH) calculation
-----------------------------------------

The self-consistent phonon mode (``MODE = SCPH``) computes temperature-dependent phonon frequencies by solving the following equation self-consistently [Tadano2015]_:

.. math::
    :label: scph_v_iter

    V_{\boldsymbol{q}ij}^{[n]} = \omega_{\boldsymbol{q}i}^{2}\delta_{ij}+\frac{1}{2}\sum_{\boldsymbol{q}_{1},k,\ell}F_{\boldsymbol{q}\boldsymbol{q}_{1},ijk\ell}\mathcal{K}_{\boldsymbol{q}_{1},k\ell}^{[n-1]}.

Here, :math:`\omega_{\boldsymbol{q}j}` is the harmonic phonon frequency and :math:`F_{\boldsymbol{q}\boldsymbol{q}_{1},ijk\ell} = \Phi(-\boldsymbol{q}i;\boldsymbol{q}j;\boldsymbol{q}_{1}k;-\boldsymbol{q}_{1}\ell)` is the reciprocal representation of fourth-order force constants computed using the harmonic eigenvectors. The updated phonon frequency in the :math:`n`\ th iteration is obtained by diagonalizing the matrix :math:`V_{\boldsymbol{q}ij}^{[n]}` as 

.. math::

    \Lambda^{[n]}_{\boldsymbol{q}} = C^{[n]\dagger}_{\boldsymbol{q}}V^{[n]}_{\boldsymbol{q}}C^{[n]}_{\boldsymbol{q}},

where :math:`\omega_{\boldsymbol{q}j}^{[n]} = (\Lambda^{[n]}_{\boldsymbol{q}jj})^{\frac{1}{2}}` and :math:`C^{[n]}_{\boldsymbol{q}}` is the unitary matrix that transforms the harmonic phonon eigenvectors into anharmonic ones as :math:`W^{[n]}_{\boldsymbol{q}} = W_{\boldsymbol{q}}C^{[n]}_{\boldsymbol{q}}`. The matrix :math:`\mathcal{K}` in Eq. :eq:`scph_v_iter` is defined as

.. math::

    \mathcal{K}_{\boldsymbol{q},ij}^{[n]} &= \alpha K_{\boldsymbol{q},ij}^{[n]} + (1-\alpha) K_{\boldsymbol{q},ij}^{[n-1]},  \\
    K_{\boldsymbol{q},ij}^{[n]} 
    &= \sum_{k} C_{\boldsymbol{q},ik}^{[n]} C_{\boldsymbol{q},jk}^{[n]*} \frac{\hbar\big[1+2n(\omega_{\boldsymbol{q}k}^{[n]})\big]}{2\omega_{\boldsymbol{q}k}^{[n]}}
    = (C^{[n]}_{\boldsymbol{q}} Q^{[n]}_{\boldsymbol{q}} C^{[n]\dagger}_{\boldsymbol{q}})_{ij},\\
    Q_{\boldsymbol{q},ij}^{[n]} 
    &= \frac{\hbar\big[1+2n(\omega_{\boldsymbol{q}i}^{[n]})\big]}{2\omega_{\boldsymbol{q}i}^{[n]}}\delta_{ij}.

:math:`\alpha` is the mixing parameter, which can be changed via the ``MIXALPHA`` tag, and :math:`n(\omega)` is the Bose-Einstein distribution function.

When ``CLASSICAL = 1`` is given in input, the expectation value of the mean square displacement is computed using classical statistics, where the :math:`Q_{\boldsymbol{q}}` matrix is replaced by

.. math::

    Q_{\boldsymbol{q},ij}^{[n]} = \frac{k_{B}T}{(\omega_{\boldsymbol{q}i}^{[n]})^{2}} \delta_{ij}.

The SCPH equation is solved on the irreducible :math:`\boldsymbol{q}` grid defined by the ``KMESH_INTERPOLATE`` tag.
The :math:`\boldsymbol{q}_{1}` grid in Eq. :eq:`scph_v_iter`, given  by the ``KMESH_SCPH`` tag,  
can be finer than the :math:`\boldsymbol{q}` grid. After the SCPH iteration converges, the code computes the anharmonic correction to the harmonic force constant :math:`\Delta D(\boldsymbol{r}(\ell))` as follows:

.. math::
    
    &\Delta D(\boldsymbol{r}(\ell)) = \frac{1}{N_{q}}\sum_{\boldsymbol{q}} \Delta D(\boldsymbol{q}) e^{-i\boldsymbol{q}\cdot\boldsymbol{r}(\ell)}, \\
    &\Delta D(\boldsymbol{q}) = D_{\mathrm{SCPH}}(\boldsymbol{q}) - D_{\mathrm{Harmonic}}(\boldsymbol{q}), \\
    &D_{\mathrm{SCPH}}(\boldsymbol{q}) = W_{\boldsymbol{q}}C_{\boldsymbol{q}}^{[n]}\Lambda_{\boldsymbol{q}}^{[n]}C_{\boldsymbol{q}}^{[n]\dagger}W_{\boldsymbol{q}}^{\dagger}.

:math:`\Delta D(\boldsymbol{r}(\ell))` is saved in ``PREFIX.scph_dfc2``.

The most computationally expensive part is the calculation of matrix elements of :math:`F_{\boldsymbol{q}\boldsymbol{q}_{1},ijk\ell}`.
When ``SELF_OFFDIAG = 0`` (default), the code only computes the elements of :math:`F_{\boldsymbol{q}\boldsymbol{q}_{1},iikk}`. 
Therefore, the computational complexity is :math:`\mathcal{O}(N_{q}^{\mathrm{irred.}}N_{q_{1}}m^{2})`.
When ``SELF_OFFDIAG = 1``, the off-diagonal elements are also calculated, and the computational complexity is :math:`\mathcal{O}(N_{q}^{\mathrm{irred.}}N_{q_{1}}m^{4})`. 

.. _formalism_relax:

Structural optimization at finite temperatures (``RELAX_STR``)
--------------------------------------------------------------

With ``MODE = SCPH`` or ``MODE = QHA`` and ``RELAX_STR > 0``, *anphon* determines the crystal structure at each
temperature by minimizing the free energy with respect to the atomic positions and, optionally, the shape of the unit
cell [Oba2019]_ [Masuki2022]_ [Masuki2023]_. The force constants are those of a fixed reference structure (the structure of the ``&cell``
field, at which the IFCs were fitted). The structure at temperature :math:`T` is described relative to it, and the IFCs
at that structure are obtained by Taylor-expanding the potential energy surface around the reference. No new
DFT calculation is needed during the optimization. The tutorials
:doc:`../tutorial_pages/bto_scph_relax` and :doc:`../tutorial_pages/zno_qha_relax` show the workflow.

Structural variables
~~~~~~~~~~~~~~~~~~~~

Two sets of variables describe the deviation from the reference structure:

* the **atomic displacements** :math:`\boldsymbol{u}_{0}(\kappa)`, uniform over all unit cells (a :math:`\Gamma`-point distortion),
  which are expanded in the harmonic normal coordinates of the reference structure at :math:`\Gamma`,

  .. math::

      u_{0,\mu}(\kappa) = \frac{1}{\sqrt{M_{\kappa}}}\sum_{j} e_{\mu}(\kappa;\boldsymbol{0}j)\, q_{0,j}.

  The three acoustic coordinates are a rigid translation and are excluded, so :math:`3N_{\kappa}-3` coordinates
  :math:`q_{0,j}` remain (written to ``PREFIX.atom_disp`` as :math:`\boldsymbol{u}_{0}`);

* the **displacement gradient tensor** :math:`u_{\mu\nu}` (``RELAX_STR = 2``), which maps a position
  :math:`\boldsymbol{r}` of the reference cell to :math:`(I+u)\boldsymbol{r}`. It is kept symmetric (the rigid rotations
  are excluded), and its six independent components are written to ``PREFIX.umn_tensor``. The energy is a function
  of the Green--Lagrange strain :math:`\eta = \mathrm{sym}(u)+\frac{1}{2}uu^{\mathrm{T}}`.

Potential energy surface at the deformed structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The potential energy is built in two steps. First, the strain dependence. The static energy of the homogeneously
strained cell is expanded with the elastic constants,

.. math::
    :label: relax_u0_strain

    U_{0}(u) = U_{\mathrm{ref}} + V\sum_{\mu\nu}C_{\mu\nu}\eta_{\mu\nu}
    + \frac{V}{2}\sum_{\mu\nu,\mu'\nu'}C_{\mu\nu,\mu'\nu'}\eta_{\mu\nu}\eta_{\mu'\nu'}
    + \frac{V}{6}\sum_{\mu\nu,\mu'\nu',\mu''\nu''}C_{\mu\nu,\mu'\nu',\mu''\nu''}\eta_{\mu\nu}\eta_{\mu'\nu'}\eta_{\mu''\nu''}
    + pV\det(I+u),

where :math:`V` is the volume of the reference cell, :math:`C_{\mu\nu}` is the residual stress of the reference structure,
:math:`C_{\mu\nu,\mu'\nu'}` and :math:`C_{\mu\nu,\mu'\nu',\mu''\nu''}` are the clamped-ion second- and third-order
(Brugger) elastic constants, and :math:`p` is the hydrostatic pressure ``STAT_PRESSURE``. Clamped-ion constants are used
because the internal coordinates are optimized explicitly. The IFCs of the strained cell are expanded in the same way,

.. math::
    :nowrap:

    \begin{align}
    \Phi_{\mu}(\kappa; u) &= \Phi_{\mu}(\kappa) + \sum_{\mu_1\nu_1}\frac{\partial \Phi_{\mu}(\kappa)}{\partial u_{\mu_1\nu_1}}u_{\mu_1\nu_1}
      + \frac{1}{2}\sum_{\mu_1\nu_1,\mu_2\nu_2}\frac{\partial^2 \Phi_{\mu}(\kappa)}{\partial u_{\mu_1\nu_1}\partial u_{\mu_2\nu_2}}u_{\mu_1\nu_1}u_{\mu_2\nu_2} + \cdots, \notag \\
    \Phi_{\mu_1\mu_2}(u) &= \Phi_{\mu_1\mu_2} + \sum_{\mu\nu}\frac{\partial \Phi_{\mu_1\mu_2}}{\partial u_{\mu\nu}}u_{\mu\nu}
      + \frac{1}{2}\sum_{\mu\nu,\mu'\nu'}\frac{\partial^2 \Phi_{\mu_1\mu_2}}{\partial u_{\mu\nu}\partial u_{\mu'\nu'}}u_{\mu\nu}u_{\mu'\nu'}, \notag \\
    \Phi_{\mu_1\mu_2\mu_3}(u) &= \Phi_{\mu_1\mu_2\mu_3} + \sum_{\mu\nu}\frac{\partial \Phi_{\mu_1\mu_2\mu_3}}{\partial u_{\mu\nu}}u_{\mu\nu}, \notag
    \end{align}

where the atomic indices are omitted for brevity, and the quartic IFCs are kept fixed (unless ``STRAIN_FC5 = 1``).
The first-order term of :math:`\Phi_{\mu}(\kappa;u)` is the strain--force coupling, the force acting on the atoms of a
uniformly strained cell. It vanishes for atoms on centrosymmetric sites, but not in general: in diamond, for instance, a
shear strain displaces the two sublattices against each other. A strain derivative of the
:math:`n`\ th-order IFCs can be computed from the :math:`(n+1)`\ th-order ones,

.. math::

    \frac{\partial \Phi_{\mu_1\cdots\mu_n}(\ell_1\kappa_1;\dots;\ell_n\kappa_n)}{\partial u_{\mu\nu}}
    = \sum_{\ell\kappa}\Phi_{\mu_1\cdots\mu_n\mu}(\ell_1\kappa_1;\dots;\ell_n\kappa_n;\ell\kappa)\, r_{\nu}(\ell\kappa),

as for the :ref:`Grüneisen parameters <formalism_gruneisen>`, or from DFT calculations of strained supercells.
Which source is used for :math:`C`, for the strain--force coupling, and for :math:`\partial\Phi_{\mu_1\mu_2}/\partial u_{\mu\nu}`
is selected by :ref:`STRAIN_COUPLING <anphon_strain_coupling>`. The terms :math:`\partial^{2}\Phi_{\mu_1\mu_2}/\partial u^{2}` and
:math:`\partial\Phi_{\mu_1\mu_2\mu_3}/\partial u` are always computed from the quartic IFCs, whereas the second and third strain
derivatives of :math:`\Phi_{\mu}(\kappa)` are included (computed from the cubic and quartic IFCs) only when the tens digit of
``STRAIN_COUPLING`` is 2, and are zero otherwise.

Second, the potential is re-expanded around the displaced atomic positions. Writing :math:`\Phi_{n}` for the strained
IFCs of order :math:`n` and :math:`\Phi_{n}\cdot\boldsymbol{u}_{0}^{k}` for the contraction of :math:`k` of their indices
with :math:`\boldsymbol{u}_{0}` (summed over all atoms of the supercell), the IFCs at the structure :math:`(u,\boldsymbol{u}_{0})` are

.. math::
    :nowrap:

    \begin{align}
    \tilde{\Phi}_{1} &= \Phi_{1} + \Phi_{2}\cdot\boldsymbol{u}_{0} + \frac{1}{2}\Phi_{3}\cdot\boldsymbol{u}_{0}^{2} + \frac{1}{6}\Phi_{4}\cdot\boldsymbol{u}_{0}^{3}, &
    \tilde{\Phi}_{2} &= \Phi_{2} + \Phi_{3}\cdot\boldsymbol{u}_{0} + \frac{1}{2}\Phi_{4}\cdot\boldsymbol{u}_{0}^{2}, \notag \\
    \tilde{\Phi}_{3} &= \Phi_{3} + \Phi_{4}\cdot\boldsymbol{u}_{0}, &
    \tilde{\Phi}_{4} &= \Phi_{4}, \notag
    \end{align}

and the static energy becomes
:math:`U_{0}(u,\boldsymbol{u}_{0}) = U_{0}(u) + \Phi_{1}\cdot\boldsymbol{u}_{0} + \frac{1}{2}\Phi_{2}\cdot\boldsymbol{u}_{0}^{2} + \frac{1}{6}\Phi_{3}\cdot\boldsymbol{u}_{0}^{3} + \frac{1}{24}\Phi_{4}\cdot\boldsymbol{u}_{0}^{4}`.
Because the potential is truncated at the fourth order in the displacements, this re-expansion is exact within the model.

Free energy and its gradients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The free energy minimized at each temperature is

.. math::

    F(u, \boldsymbol{u}_{0}, T) = U_{0}(u,\boldsymbol{u}_{0}) + F_{\mathrm{vib}}[\tilde{\Phi}_{2}, \tilde{\Phi}_{4}; T].

In ``MODE = SCPH``, :math:`F_{\mathrm{vib}}` is the SCP free energy obtained by solving the
:ref:`SCPH equation <formalism_SCPH>` with :math:`\tilde{\Phi}_{2}` in place of the harmonic IFCs. In ``MODE = QHA``, it is the
harmonic free energy :math:`N_{q}^{-1}\sum_{\boldsymbol{q}j}[\hbar\omega_{\boldsymbol{q}j}/2 + k_{B}T\ln(1-e^{-\hbar\omega_{\boldsymbol{q}j}/k_{B}T})]` of the
frequencies obtained from :math:`\tilde{\Phi}_{2}` alone. In this case the quartic IFCs enter only through the
renormalization of :math:`\tilde{\Phi}_{2}` by :math:`\boldsymbol{u}_{0}` and :math:`u`.

The SCP free energy is variational, so its derivative with respect to any parameter :math:`\lambda` of the potential is the
thermal average :math:`\langle\partial U/\partial\lambda\rangle` with the SCP (Gaussian) density matrix. The same holds
for the QHA free energy, for which :math:`\partial F_{\mathrm{vib}}/\partial\lambda = \frac{1}{2}\,\partial\tilde{\Phi}_{2}/\partial\lambda\cdot\langle\boldsymbol{u}\boldsymbol{u}\rangle`.
In both cases the gradients are

.. math::
    :label: relax_gradients

    \frac{\partial F}{\partial u_{0,\mu}(\kappa)} &= \tilde{\Phi}_{\mu}(\kappa)
    + \frac{1}{2}\sum_{\ell_1\kappa_1\nu,\ell_2\kappa_2\lambda}\tilde{\Phi}_{\mu\nu\lambda}(0\kappa;\ell_1\kappa_1;\ell_2\kappa_2)
    \langle u_{\nu}(\ell_1\kappa_1)u_{\lambda}(\ell_2\kappa_2)\rangle, \\
    \frac{\partial F}{\partial u_{\mu\nu}} &= \frac{\partial U_{0}(u,\boldsymbol{u}_{0})}{\partial u_{\mu\nu}}
    + \frac{1}{2}\sum_{\ell_1\kappa_1\mu_1,\ell_2\kappa_2\mu_2}\frac{\partial\tilde{\Phi}_{\mu_1\mu_2}(\ell_1\kappa_1;\ell_2\kappa_2)}{\partial u_{\mu\nu}}
    \langle u_{\mu_1}(\ell_1\kappa_1)u_{\mu_2}(\ell_2\kappa_2)\rangle,

i.e. the force is the static force plus the thermal average of the cubic term, and the stress is the static stress
plus the vibrational (thermal and zero-point) stress. With ``STRAIN_FC5 = 1`` the stress acquires the additional term
:math:`\frac{1}{8}\,\partial\Phi_{4}/\partial u_{\mu\nu}\cdot\langle\boldsymbol{u}\boldsymbol{u}\rangle\langle\boldsymbol{u}\boldsymbol{u}\rangle`.
The displacement correlation function is evaluated from the SCP (or QHA) phonons on the ``KMESH_SCPH``
(``KMESH_QHA``) mesh,

.. math::

    \langle u_{\mu}(\ell\kappa)u_{\nu}(\ell'\kappa')\rangle = \frac{1}{N_{q}\sqrt{M_{\kappa}M_{\kappa'}}}\sum_{\boldsymbol{q}j}
    \frac{\hbar\big[1+2n(\omega_{\boldsymbol{q}j})\big]}{2\omega_{\boldsymbol{q}j}}\,
    e_{\mu}(\kappa;\boldsymbol{q}j)\,e^{*}_{\nu}(\kappa';\boldsymbol{q}j)\,
    e^{\mathrm{i}\boldsymbol{q}\cdot(\boldsymbol{r}(\ell)-\boldsymbol{r}(\ell'))},

which is the matrix :math:`K_{\boldsymbol{q}}` of the SCPH section expressed in the atomic basis. When ``CLASSICAL = 1``
the factor :math:`\hbar(1+2n)/2\omega` is replaced by :math:`k_{B}T/\omega^{2}`. Imaginary modes on this mesh are
excluded from the free energy but enter the gradients with :math:`|\omega|`, and a warning is printed.

The equilibrium structure satisfies :math:`\partial F/\partial q_{0,j} = 0` for all optical coordinates and, with
``RELAX_STR = 2``, :math:`\partial F/\partial u_{\mu\nu} = 0`. A symmetric structure is a stationary point by symmetry,
so a spontaneous symmetry breaking (e.g. the ferroelectric transition of BaTiO\ :sub:`3`) can only be obtained from an
initial structure that is already slightly distorted (``&displace`` field, ``SET_INIT_STR``).

Optimization algorithm
~~~~~~~~~~~~~~~~~~~~~~

At each temperature, the following steps are repeated until the changes of :math:`\boldsymbol{u}_{0}` and :math:`u` fall below
``COORD_CONV_TOL`` and ``CELL_CONV_TOL`` (and, when given, the gradients below ``GRADIENT_CONV_TOL`` and
``CELL_GRADIENT_CONV_TOL``), or ``MAX_STR_ITER`` is reached:

1. the IFCs :math:`\tilde{\Phi}_{n}` are evaluated at the current structure;
2. the SCPH equation is solved (``MODE = SCPH``), or the dynamical matrix of :math:`\tilde{\Phi}_{2}` is diagonalized (``MODE = QHA``);
3. the gradients :eq:`relax_gradients` are computed;
4. the structure is updated.

The update uses the state vector :math:`\boldsymbol{x} = (q_{0,j}, u_{\mu\nu})` and the gradient
:math:`\boldsymbol{g} = \partial F/\partial\boldsymbol{x}`. With the Newton-like method (``RELAX_ALGO = 2``, default),
:math:`\boldsymbol{x}\leftarrow\boldsymbol{x}-\beta H^{-1}\boldsymbol{g}`. Here the coordinate block of :math:`H` is the
:math:`\Gamma`-point dynamical matrix of the current SCP (or QHA) phonons in the basis of :math:`q_{0,j}`, plus
:math:`\omega_{\mathrm{add}}^{2}` on the diagonal (``ADD_HESS_DIAG``) to keep it positive near a soft-mode instability, and the
cell block consists of the second-order elastic constants :math:`V C_{\mu\nu,\mu'\nu'}`. The coupling between the two blocks
is neglected, unless ``BUBBLE_HESS = 1`` provides the full curvature. The mixing factors :math:`\beta` are ``MIXBETA_COORD`` and
``MIXBETA_CELL``. ``RELAX_ALGO = 3`` starts a BFGS update from the same :math:`H` and accelerates it with GDIIS, and
``RELAX_ALGO = 1`` is a plain steepest descent of the internal coordinates only; it keeps the cell fixed even when
``RELAX_STR = 2``.

Variants in ``MODE = QHA``
~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``QHA_SCHEME`` tag selects how the conditions above are imposed when the cell is relaxed (``RELAX_STR = 2``) [Masuki2023]_:

* ``QHA_SCHEME = 0`` minimizes :math:`F` with respect to both :math:`\boldsymbol{u}_{0}` and :math:`u` (full QHA).
* ``QHA_SCHEME = 1`` (zero static internal stress approximation, ZSISA): the internal coordinates are determined by the
  static condition :math:`\partial U_{0}/\partial\boldsymbol{q}_{0} = 0` at every strain, and only the cell minimizes the free energy.
  The cell gradient is the total derivative
  :math:`\mathrm{d}F/\mathrm{d}u_{\mu\nu} = \partial F/\partial u_{\mu\nu} + \sum_{j}(\partial F/\partial q_{0,j})(\partial q_{0,j}/\partial u_{\mu\nu})`,
  with the internal-strain response :math:`\partial\boldsymbol{q}_{0}/\partial u_{\mu\nu} = -\tilde{\Phi}_{2}(\Gamma)^{-1}\,\partial\tilde{\Phi}_{1}/\partial u_{\mu\nu}`.
* ``QHA_SCHEME = 2`` (volumetric ZSISA, v-ZSISA): only the volume is determined by the free energy. The cell shape and the
  internal coordinates at a given volume follow from the static energy, i.e. the vibrational stress acts only along
  the strain direction that changes :math:`\det(I+u)` at the lowest static-energy cost (obtained with the ZSISA elastic
  constants).

``RELAX_STR = 3`` (``MODE = QHA`` only) is the lowest-order perturbative solution of the same problem: the gradients of
:math:`F_{\mathrm{vib}}` are evaluated once with the harmonic phonons of the reference structure, and the structural change is
obtained from a single linear solve with the static Hessian (optical :math:`\Gamma` modes and elastic constants), without
iteration.

``RELAX_STR = 4`` keeps the strain fixed at the value of the ``&strain`` field and relaxes only :math:`\boldsymbol{u}_{0}`. The
resulting :math:`\bar{F}(u,T) = \min_{\boldsymbol{u}_{0}}F(u,\boldsymbol{u}_{0},T)` can be differentiated numerically with respect
to the strain to obtain the isothermal elastic constants at finite temperature.

The relaxed structure (and the IFCs :math:`\tilde{\Phi}_{3}`, :math:`\tilde{\Phi}_{4}` at that structure) can be used in
subsequent phonon and thermal conductivity calculations through
:ref:`RELAXED_STRUCTURE <anphon_relaxed_structure>`.


````

.. [Parlinski1998] K\. Parlinski, Z. Q. Li, and Y. Kawazoe, Phys. Rev. Lett. **81**, 3298 (1998).

.. [Wang2010] Y\. Wang *et al.*, J. Phys.: Condens. Matter **22**, 202201 (2010).

.. [Hafner1993] J\. Hafner and M. Krajci, J. Phys.: Condens. Matter **5**, 2489 (1993).

.. [Tamura1983] S\. -I. Tamura, Phys. Rev. B **27**, 858 (1983).

.. [Blochl1994] P\. E. Bl\ |umulaut_o|\ chl, O. Jepsen, and O. K. Andersen, Phys. Rev. B **49**, 16223 (1994).

.. [Simoncelli2019] M\. Simoncelli, N. Marzari, and F. Mauri, Nat. Phys. **15**, 809 (2019).

.. [Allen1993] P\. B. Allen and J. L. Feldman, Phys. Rev. B **48**, 12581 (1993).

.. [Gonze1997] X\. Gonze and C. Lee, Phys. Rev. B **55**, 10355 (1997).

.. [Simoncelli2022] M\. Simoncelli, N. Marzari, and F. Mauri, Phys. Rev. X **12**, 041011 (2022).

.. [Omini1996] M\. Omini and A. Sparavigna, Phys. Rev. B **53**, 9064 (1996).

.. [Fugallo2013] G\. Fugallo, M. Lazzeri, L. Paulatto, and F. Mauri, Phys. Rev. B **88**, 045430 (2013).

.. [Chaput2013] L\. Chaput, Phys. Rev. Lett. **110**, 265506 (2013).

.. [Yates2007] J\. R. Yates, X. Wang, D. Vanderbilt, and I. Souza, Phys. Rev. B **75**, 195121 (2007).

.. [Feng2016] T\. Feng and X. Ruan, Phys. Rev. B **93**, 045202 (2016).
