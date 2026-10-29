import os
import shutil

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("ase")
pytest.importorskip("spglib")
pytest.importorskip("lxml")

from fcsio.reader import ForceConstantParser  # noqa: E402

SI_H5 = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "example",
    "Si",
    "anharm_IFCs",
    "2_generate_config",
    "si222_harmonic.h5",
)


def test_h5_without_primcell(tmp_path):
    """alm run without PRIMCELL stores the supercell as /PrimitiveCell."""
    fname = str(tmp_path / "si222_noprim.h5")
    shutil.copy(SI_H5, fname)
    with h5py.File(fname, "r+") as f:
        del f["PrimitiveCell"]
        f.copy("SuperCell", "PrimitiveCell")

    ref = ForceConstantParser(SI_H5)
    parser = ForceConstantParser(fname)
    assert len(parser.primitive_cell[1]) == 2
    np.testing.assert_array_equal(
        parser.transformation_matrix_int, ref.transformation_matrix_int
    )
    fc_ref = ref.read(maxorder=2)["fc2"]
    fc = parser.read(maxorder=2)["fc2"]
    for a, b in zip(fc, fc_ref):
        np.testing.assert_allclose(a, b)
