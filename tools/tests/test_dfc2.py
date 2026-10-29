"""dfc2.py with a PREFIX.scph.h5 state file: the stored renormalized FC2 is written as is."""

import h5py
import numpy as np

import dfc2


def test_state_file_written_directly(scph_state_h5, tmp_path):
    path, _, base, tdep = scph_state_h5
    out = tmp_path / "fc2_200K.h5"
    dfc2.write_from_state_file(str(path), str(out), 200.0)
    with h5py.File(out, "r") as f:
        np.testing.assert_array_equal(
            f["ForceConstants/Order2/force_constant_values"][:], tdep[1]
        )
        assert (
            f["ForceConstants/Order2/force_constant_values"].attrs["unit"]
            == "Ry/bohr^2"
        )
        assert "ForceConstants/Order2_temperature_dependent" not in f
        assert "settings" not in f and "mode" not in f.attrs
        assert f.attrs["schema"] == "alamode:force_constants"
        assert f["SuperCell/number_of_atoms"][()] == 4
