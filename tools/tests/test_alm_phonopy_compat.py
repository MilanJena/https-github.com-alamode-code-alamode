"""The in-repo alm package must work with phonopy's ALM interface unchanged.

``_phonopy_sequence`` replays the calls of phonopy 4.7.2
(phonopy/interface/alm.py: ``_update_options``, ``ALMFCSolver._run`` and
``_extract_fc_from_alm``) without importing phonopy.  The model is 8-atom
diamond Si with isotropic nearest-neighbour springs, so the exact FC2 is known.
"""

import ctypes
import os

import numpy as np
import pytest

A_LAT = 10.26  # bohr
K_SPRING = 0.1  # Ry/bohr^2
POSITIONS = np.array(
    [
        [0.00, 0.00, 0.00],
        [0.00, 0.50, 0.50],
        [0.50, 0.00, 0.50],
        [0.50, 0.50, 0.00],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ]
)
NUMBERS = [14] * 8


def _exact_fc2():
    nat = len(POSITIONS)
    fc2 = np.zeros((nat, nat, 3, 3))
    for i in range(nat):
        for j in range(nat):
            d = POSITIONS[j] - POSITIONS[i]
            d = (d - np.rint(d)) * A_LAT
            if i != j and np.isclose(np.linalg.norm(d), A_LAT * np.sqrt(3) / 4):
                e = d / np.linalg.norm(d)
                fc2[i, j] -= K_SPRING * np.outer(e, e)
                fc2[i, i] += K_SPRING * np.outer(e, e)
    return fc2


def _dataset(nsnap=12, amp=0.02):
    fc2 = _exact_fc2()
    nat = len(POSITIONS)
    u = np.random.default_rng(0).normal(scale=amp, size=(nsnap, nat, 3))
    f = -np.einsum("ijab,njb->nia", fc2, u)
    return u, f


def _phonopy_sequence(maxorder, cutoff_radii, disps, forces, options):
    """phonopy 4.7.2's ALM calls, in order (log_level=0, iconst=0)."""
    from alm import ALM, optimizer_control_data_types

    # _update_options: coerce "key = value" strings with the type map
    alm_option_types = dict(optimizer_control_data_types)
    alm_options = {
        "solver": "dense",
        "nbody": None,
        "symmetrization_basis": "Lattice",
        "output_filename_prefix": None,
    }
    for key, val in options.items():
        alm_options[key] = alm_option_types[key](val)

    lattice = np.eye(3) * A_LAT
    natom = len(NUMBERS)
    with ALM(lattice, POSITIONS, NUMBERS) as alm:
        _ = alm.kind_names.values()
        alm.output_filename_prefix = alm_options["output_filename_prefix"]
        alm.verbosity = 0
        alm.define(
            maxorder,
            cutoff_radii=cutoff_radii,
            nbody=alm_options["nbody"],
            symmetrization_basis=alm_options["symmetrization_basis"],
        )
        alm.displacements = disps
        alm.forces = forces
        optcontrol = {}
        for key in optimizer_control_data_types:
            if key in alm_options:
                optcontrol[key] = alm_options[key]
        if optcontrol:
            alm.optimizer_control = optcontrol
            if "cross_validation" in optcontrol and optcontrol["cross_validation"] > 0:
                alm.optimize(solver=alm_options["solver"])
                optcontrol["cross_validation"] = 0
                optcontrol["l1_alpha"] = alm.cv_l1_alpha
                alm.optimizer_control = optcontrol
        alm.set_constraint(translation=False)  # iconst = 0
        alm.optimize(solver=alm_options["solver"])

        # _extract_fc_from_alm: compact (origin) and full (all) forms
        p2s_map = alm.getmap_primitive_to_supercell()[0]
        fcs = []
        for order in range(1, maxorder + 1):
            shape = (natom,) * (order + 1) + (3,) * (order + 1)
            fc_origin = np.zeros((len(p2s_map),) + shape[1:])
            for fc_elem, indices in zip(
                *alm.get_fc(order, mode="origin", permutation=1), strict=True
            ):
                v = indices // 3
                c = indices % 3
                p = int(np.where(p2s_map == v[0])[0][0])
                fc_origin[(p,) + tuple(v[1:]) + tuple(c)] = fc_elem
            fc = np.zeros(shape)
            for fc_elem, indices in zip(
                *alm.get_fc(order, mode="all", permutation=1), strict=True
            ):
                v = indices // 3
                c = indices % 3
                fc[tuple(v) + tuple(c)] = fc_elem
            assert np.allclose(fc_origin, fc[p2s_map], atol=1e-12)
            fcs.append(fc)
    return fcs


def _flushed_output(capfd):
    ctypes.CDLL(None).fflush(None)  # C stdio buffer of the C++ library
    captured = capfd.readouterr()
    return captured.out + captured.err


def test_optimizer_control_data_types(alm_mod):
    from alm import ALM, optimizer_control_data_types

    assert "optimizer_control_data_types" in alm_mod.__all__
    assert len(optimizer_control_data_types) == 15
    with ALM(np.eye(3) * A_LAT, POSITIONS, NUMBERS) as alm:
        alm.define(1)
        oc = alm.optimizer_control
    for key, typ in optimizer_control_data_types.items():
        assert type(oc[key]) is typ, key


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize("maxorder", [1, 2])
def test_phonopy_call_sequence(alm_mod, capfd, maxorder):
    _flushed_output(capfd)  # drop output buffered by earlier tests
    u, f = _dataset()
    cutoff_radii = -np.ones((maxorder, 1, 1))
    if maxorder == 2:
        cutoff_radii[1] = 5.0  # FC3: nearest neighbours only
    fcs = _phonopy_sequence(
        maxorder,
        cutoff_radii,
        u,
        f,
        {"linear_model": "1", "tolerance_iteration": "1e-8"},
    )
    assert _flushed_output(capfd) == ""
    assert np.abs(fcs[0] - _exact_fc2()).max() < 1e-10
    if maxorder == 2:
        assert np.abs(fcs[1]).max() < 1e-10


@pytest.mark.filterwarnings("error")
def test_phonopy_run_alm(alm_mod, capfd):
    """Through phonopy itself. Optional locally; CI sets ALM_REQUIRE_PHONOPY=1
    (no skip: CI fails on skips)."""
    _flushed_output(capfd)  # drop output buffered by earlier tests
    try:
        import phonopy  # noqa: F401
    except ModuleNotFoundError:
        if os.environ.get("ALM_REQUIRE_PHONOPY") == "1":
            raise
        return
    from phonopy.interface.alm import run_alm
    from phonopy.structure.atoms import PhonopyAtoms
    from phonopy.structure.cells import Primitive

    cell = PhonopyAtoms(
        numbers=NUMBERS, cell=np.eye(3) * A_LAT, scaled_positions=POSITIONS
    )
    primitive = Primitive(cell, np.eye(3))
    u, f = _dataset()
    for options in (None, "iconst = 0, linear_model = 1"):
        fcs = run_alm(cell, primitive, u, f, 1, options=options)
        assert _flushed_output(capfd) == ""
        assert np.abs(fcs[2] - _exact_fc2()).max() < 1e-10
