"""Checks of the anphon Python layer (python/anphon/_api.py): run() from dicts vs
the anphon CLI (Si and BaTiO3 bands, Si kappa.h5), Phonon.frequencies vs the CLI band
structure, Result.h5("kappa"), errors raised as RuntimeError with the process still usable
afterwards, and mpirun -np 2 smoke runs (test names containing "mpi").

All serial anphon runs from Python happen in one child process, so that a fatal error
cannot take down the test session; the parent only compares files."""

import json
import os
import shutil
import signal
import subprocess
import sys

import numpy as np
import pytest

pytest.importorskip("anphon")
pytest.importorskip("h5py")

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ANPHON_BIN = os.path.join(REPO, "_build", "anphon", "anphon")
SI = os.path.join(REPO, "example", "Si")
SI_XML = os.path.join(
    SI, "anharm_IFCs", "1_harmonic", "reference", "si222_harmonic.xml"
)
SI_CUBIC = os.path.join(SI, "reference", "si222_cubic.h5")
BTO = os.path.join(REPO, "example", "BaTiO3")
BTO_XML = os.path.join(BTO, "anharm_IFCs", "cBTO222_harmonic.xml")
BTO_BORN = os.path.join(BTO, "scph_relax", "BORNINFO")
# attributes that identify the build or the moment of the run, not the result
KAPPA_VOLATILE = {"@created_date", "@alamode_git_commit", "@alamode_version"}

A_SI = 10.2631026
FCC = [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
PATH = [
    ("G", [0, 0, 0], "X", [0.5, 0.5, 0], 51),
    ("X", [0.5, 0.5, 1], "G", [0, 0, 0], 51),
    ("G", [0, 0, 0], "L", [0.5, 0.5, 0.5], 51),
]
# rows of the .bands file at Gamma, X = (1/2, 1/2, 0) and L = (1/2, 1/2, 1/2)
QPOINTS = {0: [0, 0, 0], 50: [0.5, 0.5, 0], 152: [0.5, 0.5, 0.5]}

CLI_INPUT = f"""&general
  PREFIX = si
  MODE = phonons
  FCSFILE = {SI_XML}
/
&cell
  {A_SI}
  0.0 0.5 0.5
  0.5 0.0 0.5
  0.5 0.5 0.0
/
&kpoint
  1
  G 0 0 0 X 0.5 0.5 0 51
  X 0.5 0.5 1 G 0 0 0 51
  G 0 0 0 L 0.5 0.5 0.5 51
/
"""

_CHILD = """
import json, os, sys
import numpy as np
import anphon

out, cfg = sys.argv[1], json.loads(sys.argv[2])
cell = dict(a=cfg["a"], lattice=cfg["fcc"])

# (i) run() from dicts
r = anphon.run(
    general=dict(PREFIX="si", MODE="phonons", FCSFILE=cfg["xml"]),
    cell=cell,
    kpoint=dict(KPMODE=1, path=cfg["path"]),
    workdir=os.path.join(out, "py"),
)
assert r.prefix == "si" and "si.bands" in r.files(), r.files()
try:
    r.h5("kappa")
except FileNotFoundError as e:
    assert "si.kappa.h5" in str(e)
else:
    raise AssertionError("no FileNotFoundError")

# run() on an existing input file
r2 = anphon.run(os.path.join(out, "cli", "si.in"), workdir=os.path.join(out, "py2"))
assert r2.prefix == "si" and r2.path(".bands").is_file()

# (ii) Phonon.frequencies; the markers delimit its log (quiet by default)
print("<<QUIET", flush=True)
ph = anphon.Phonon(cfg["xml"], cell=cell)
w, e = ph.frequencies(cfg["q"], eigenvectors=True)
np.save(os.path.join(out, "w.npy"), w)
np.save(os.path.join(out, "e.npy"), e)
assert np.array_equal(ph.frequencies(cfg["q"]), w)
print("QUIET>>\\n<<LOUD", flush=True)
anphon.Phonon(cfg["xml"], cell=cell, general=dict(verbosity=1)).frequencies(cfg["q"][0])
print("LOUD>>", flush=True)

# (iii) a tiny kappa run
if cfg["cubic"]:
    k = anphon.run(
        general=dict(PREFIX="si", MODE="kappa", FCSFILE=cfg["cubic"], KD=["Si"],
                     MASS=28.0855, TMIN=300, TMAX=300),
        cell=dict(a=10.203, lattice=cfg["fcc"]),
        kpoint=dict(KPMODE=2, mesh=[4, 4, 4]),
        kappa=dict(SOLVER="RTA"),
        workdir=os.path.join(out, "kappa"),
    ).h5("kappa")
    assert k["@schema"] == "alamode:kappa_result", k["@schema"]
    np.save(os.path.join(out, "kappa.npy"), k["kappa"]["kappa_peierls"])

# (iv) another system (BaTiO3, non-analytic correction) in the same process; the
# next Si run must not see any of its state
anphon.run(
    general=dict(PREFIX="bto", MODE="phonons", FCSFILE=cfg["bto_xml"], KD=["Ba", "Ti", "O"],
                 BORNINFO=cfg["bto_born"], NONANALYTIC=2),
    cell=dict(a=7.53159676409, lattice=[[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    kpoint=dict(KPMODE=1, path=[("G", [0, 0, 0], "X", [0.5, 0, 0], 11),
                                ("X", [0.5, 0, 0], "M", [0.5, 0.5, 0], 11),
                                ("M", [0.5, 0.5, 0], "R", [0.5, 0.5, 0.5], 11)]),
    workdir=os.path.join(out, "bto"),
)
anphon.run(os.path.join(out, "cli", "si.in"), workdir=os.path.join(out, "py3"))

# (v) errors raise RuntimeError, leave the cwd alone, and the next run works
import h5py

bad = os.path.join(out, "bad")
os.makedirs(bad)
with open(os.path.join(bad, "garbage.xml"), "w") as f:
    f.write("this is <not xml")
with h5py.File(os.path.join(bad, "noschema.h5"), "w") as f:
    f.attrs["schema"] = "nonsense"
    f["Foo/bar"] = [1.0, 2.0]
with open(os.path.join(bad, "BORNINFO"), "w") as f:
    f.write("13.0 0 0  0 13.0 0  0 0 13.0")  # eps_inf only, no Z*
si = dict(PREFIX="bad", MODE="phonons", FCSFILE=cfg["xml"])
cases = {
    "xml": (dict(si, FCSFILE=os.path.join(bad, "garbage.xml")), "Cannot read FCSFILE"),
    # anphon does not check the schema of an FCSFILE, only that the groups it needs exist
    "h5": (dict(si, FCSFILE=os.path.join(bad, "noschema.h5")), "has no /SuperCell group"),
    "born": (dict(si, NONANALYTIC=1, BORNINFO=os.path.join(bad, "BORNINFO")),
             "has 9 numbers"),
}
errors = {}
for name, (general, fragment) in cases.items():
    cwd = os.getcwd()
    try:
        anphon.run(general=general, cell=cell, kpoint=dict(KPMODE=1, path=cfg["path"]),
                   workdir=os.path.join(bad, name))
    except RuntimeError as e:
        errors[name] = [type(e).__name__, str(e)]
        assert fragment in str(e), (name, str(e))
    else:
        raise AssertionError(f"no RuntimeError for {name}")
    assert os.getcwd() == cwd, (name, os.getcwd())
    anphon.run(os.path.join(out, "cli", "si.in"), workdir=os.path.join(out, "after_" + name))
with open(os.path.join(out, "errors.json"), "w") as f:
    json.dump(errors, f)
print("CHILD_OK")
"""

# mpirun -np 2: each rank prints OK<rank>/<size> once it sees the complete outputs
_MPI_CHILD = """
import os, sys
if sys.argv[1] == "mpi4py_first":
    from mpi4py import MPI
import anphon
if sys.argv[1] == "anphon_first":
    from mpi4py import MPI  # finds MPI initialised by anphon
from anphon._anphon import comm_size

rank = int(next(os.environ[k] for k in ("OMPI_COMM_WORLD_RANK", "PMIX_RANK", "PMI_RANK")
                if k in os.environ))
src, wd = sys.argv[2], sys.argv[3]
if sys.argv[1] == "fail" and rank == 1:
    src = os.path.join(wd, "missing.in")  # rank-local failure while preparing
r = anphon.run(src, workdir=wd)
assert r.path(".bands").is_file(), r.files()
if sys.argv[1].endswith("_first"):
    assert MPI.COMM_WORLD.Get_size() == comm_size() and MPI.COMM_WORLD.allreduce(1) == 2
    assert MPI.Is_initialized() and not MPI.Is_finalized()
print(f"OK{rank}/{comm_size()}", flush=True)
"""


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    if not os.path.exists(ANPHON_BIN):
        pytest.skip("anphon binary (_build/anphon/anphon) not found")
    out = tmp_path_factory.mktemp("anphon_api")
    # one thread: the OpenMP reductions of a kappa run differ in the last bit otherwise
    env = dict(os.environ, OMP_NUM_THREADS="1")
    cli = out / "cli"
    cli.mkdir()
    (cli / "si.in").write_text(CLI_INPUT)
    with open(cli / "si.log", "w") as log:
        subprocess.run(
            [ANPHON_BIN, "si.in"], cwd=cli, stdout=log, check=True, env=env, timeout=600
        )

    cubic = SI_CUBIC if os.path.exists(SI_CUBIC) else ""
    cfg = dict(
        a=A_SI,
        fcc=FCC,
        xml=SI_XML,
        path=PATH,
        q=list(QPOINTS.values()),
        cubic=cubic,
        bto_xml=BTO_XML,
        bto_born=BTO_BORN,
    )
    res = subprocess.run(
        [sys.executable, "-c", _CHILD, str(out), json.dumps(cfg)],
        capture_output=True,
        text=True,
        env=env,
        timeout=600,
    )
    assert res.returncode == 0 and "CHILD_OK" in res.stdout, (
        res.stdout[-4000:] + res.stderr[-4000:]
    )
    (out / "child.log").write_text(res.stdout)
    # the CLI on the input files the child's run() wrote
    for name in ("bto", "kappa") if cubic else ("bto",):
        inp = next((out / name).glob("*.in"))
        (out / f"cli_{name}").mkdir()
        shutil.copy(inp, out / f"cli_{name}")
        with open(out / f"cli_{name}" / "cli.log", "w") as log:
            subprocess.run(
                [ANPHON_BIN, inp.name],
                cwd=out / f"cli_{name}",
                stdout=log,
                check=True,
                env=env,
                timeout=600,
            )
    return out


def test_run_from_dicts_matches_cli(runs):
    cli = (runs / "cli" / "si.bands").read_text()
    for d in ("py", "py2", "py3"):  # py3: after a BaTiO3 run in the same process
        assert (runs / d / "si.bands").read_text() == cli, d


def test_other_system_matches_cli(runs):
    py = (runs / "bto" / "bto.bands").read_text()
    assert py == (runs / "cli_bto" / "bto.bands").read_text()
    assert len(np.loadtxt(runs / "bto" / "bto.bands")) == 33


def test_kappa_h5_matches_cli(runs):
    from anphon import read_h5

    if not (runs / "cli_kappa").exists():
        pytest.skip("example/Si/reference/si222_cubic.h5 not found")
    py = read_h5(runs / "kappa" / "si.kappa.h5")
    cli = read_h5(runs / "cli_kappa" / "si.kappa.h5")
    assert set(py) == set(cli) and KAPPA_VOLATILE <= set(py)

    def compare(a, b, path):
        assert set(a) == set(b), path
        for k in a:
            if path == "" and k in KAPPA_VOLATILE:
                continue
            if isinstance(a[k], dict):
                compare(a[k], b[k], f"{path}/{k}")
            elif isinstance(a[k], (str, list)):
                assert a[k] == b[k], f"{path}/{k}"
            else:
                assert np.array_equal(a[k], b[k]), f"{path}/{k}"

    compare(py, cli, "")


def test_errors_raise_and_next_run_works(runs):
    errors = json.loads((runs / "errors.json").read_text())
    assert sorted(errors) == ["born", "h5", "xml"]
    assert all(t == "RuntimeError" for t, _ in errors.values()), errors
    cli = (runs / "cli" / "si.bands").read_text()
    for name in errors:
        assert (runs / f"after_{name}" / "si.bands").read_text() == cli, name


def test_phonon_frequencies_match_bands(runs):
    bands = np.loadtxt(runs / "cli" / "si.bands")[:, 1:]
    w = np.load(runs / "w.npy")
    assert w.shape == (3, 6)
    # .bands carries 7 significant digits
    assert np.allclose(w, bands[list(QPOINTS)], rtol=2e-6, atol=1e-6)
    assert np.all(np.abs(w[0, :3]) < 1e-4)  # acoustic modes at Gamma


def test_phonon_log(runs):
    log = (runs / "child.log").read_text()
    quiet = log.split("<<QUIET\n")[1].split("QUIET>>")[0]
    loud = log.split("<<LOUD\n")[1].split("LOUD>>")[0]
    assert quiet == "", quiet  # VERBOSITY = 0 by default, no banner
    assert "Program ANPHON" in loud and "Job finished" in loud, loud
    assert "Program ANPHON" in log.split("<<QUIET")[0]  # run() keeps the CLI default


def test_phonon_eigenvectors_orthonormal(runs):
    e = np.load(runs / "e.npy")
    assert e.shape == (3, 6, 6) and np.iscomplexobj(e)
    for eq in e:  # rows are the modes
        assert np.allclose(eq @ eq.conj().T, np.eye(6), atol=1e-10)
    # axis 1 is the mode: at Gamma the acoustic modes move both atoms alike,
    # the optical ones in antiphase (fails if the mode/component axes were swapped)
    assert np.allclose(e[0, :3, :3], e[0, :3, 3:], atol=1e-8)
    assert np.allclose(e[0, 3:, :3], -e[0, 3:, 3:], atol=1e-8)


def test_result_h5_kappa(runs):
    if not (runs / "kappa.npy").exists():
        pytest.skip("example/Si/reference/si222_cubic.h5 not found")
    kap = np.load(runs / "kappa.npy")
    assert kap.shape == (1, 3, 3)
    assert np.all(np.diag(kap[0]) > 0)
    assert np.allclose(kap[0], kap[0, 0, 0] * np.eye(3), rtol=1e-6, atol=1e-8)


# mpirun -np 2 (each run in its own process group, killed on timeout)


def _mpirun():
    exe = os.environ.get("ANPHON_TEST_MPIRUN")
    if not exe:
        exe = os.path.join(os.path.dirname(sys.executable), "mpirun")
        exe = exe if os.access(exe, os.X_OK) else shutil.which("mpirun")
    if not exe:
        pytest.skip("mpirun not found (set ANPHON_TEST_MPIRUN)")
    return exe


def _run_mpi(mode, runs, tmp_path, timeout=120):
    cmd = [_mpirun(), "-np", "2", sys.executable, "-c", _MPI_CHILD, mode]
    cmd += [str(runs / "cli" / "si.in"), str(tmp_path)]
    env = dict(os.environ, OMP_NUM_THREADS="1")
    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env,
        start_new_session=True,
    )
    try:
        out, _ = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        # SIGTERM first so mpirun can take down ranks outside its process group
        for sig in (signal.SIGTERM, signal.SIGKILL):
            try:
                os.killpg(p.pid, sig)
                out, _ = p.communicate(timeout=10)
                break
            except ProcessLookupError:
                out = ""
                break
            except subprocess.TimeoutExpired:
                out = ""
        pytest.fail(f"mpirun {mode} timed out after {timeout} s:\n{out[-4000:]}")
    return p.returncode, out


@pytest.mark.parametrize("mode", ["anphon", "anphon_first", "mpi4py_first"])
def test_mpi_run(mode, runs, tmp_path):
    if mode != "anphon":
        pytest.importorskip("mpi4py")
    rc, out = _run_mpi(mode, runs, tmp_path)
    assert rc == 0 and "OK0/2" in out and "OK1/2" in out, out[-4000:]
    cli = (runs / "cli" / "si.bands").read_text()
    assert (tmp_path / "si.bands").read_text() == cli


def test_mpi_rank_local_failure_aborts(runs, tmp_path):
    rc, out = _run_mpi("fail", runs, tmp_path)
    assert rc != 0 and "OK0/2" not in out, out[-4000:]
    assert "anphon.run:" in out and "missing.in" in out, out[-4000:]


# pure-Python checks (no anphon run)


def test_phonon_rejects_own_general_tags():
    import anphon

    for key in (
        "prefix",
        "MODE",
        "FcsFile",
        "BORNINFO",
        "File_Format",
        "Nonanalytic",
        "KD",
    ):
        with pytest.raises(ValueError):
            anphon.Phonon(SI_XML, general={key: "x"})  # kd= etc. are arguments
    anphon.Phonon(SI_XML, general=dict(MASS=28.0855, TOLERANCE=1e-5))
    with pytest.raises(TypeError):  # no **kwargs sink for &general tags
        anphon.Phonon(SI_XML, VERBOSITY=0)


def test_phonon_copies_settings_and_checks_q(monkeypatch):
    import anphon
    import anphon._api as api

    general = dict(MASS=np.array([28.0855]))
    lattice = np.array(FCC, dtype=float)
    kd = ["Si"]
    cell = dict(a=A_SI, lattice=lattice)
    ph = anphon.Phonon(SI_XML, cell=cell, kd=kd, general=general)
    general["MODE"] = "kappa"  # changes after construction have no effect,
    general["MASS"][0] = 1.0  # also inside nested objects (deep copy)
    cell["a"] = 1.0
    lattice[0, 1] = 9.0
    kd[0] = "Ge"
    seen = {}

    class Stop(Exception):
        pass

    def fake_run(**sections):
        seen.update(sections)
        raise Stop  # before anphon would run

    monkeypatch.setattr(api, "run", fake_run)
    with pytest.raises(Stop):
        ph.frequencies([0, 0, 0])
    assert (
        seen["general"]["MODE"] == "phonons" and seen["general"]["MASS"][0] == 28.0855
    )
    assert seen["cell"]["a"] == A_SI
    np.testing.assert_array_equal(seen["cell"]["lattice"], FCC)
    assert seen["general"]["KD"] == ["Si"]
    for bad in (False, [("MASS", 1.0)], "MASS=1"):
        with pytest.raises(TypeError):
            anphon.Phonon(SI_XML, cell=cell, general=bad)
    with pytest.raises(ValueError, match="lacks"):
        anphon.format_input(general=dict(PREFIX="x"), cell={})
    assert seen["kpoint"]["qpoints"].shape == (1, 3)
    for bad in (np.zeros((2, 6)), np.zeros(6), np.zeros((1, 1, 3))):
        with pytest.raises(ValueError):
            ph.frequencies(bad)
    with pytest.raises(AttributeError):
        ph.fcsfile = "x"


def test_public_names():
    import anphon
    import anphon._api as api

    assert "run_input_file" not in dir(anphon)
    assert callable(anphon._anphon.run_input_file)  # private low-level binding
    assert set(anphon.__all__) == set(api.__all__)
    assert set(anphon.__all__) == {"run", "Result", "Phonon", "format_input", "read_h5"}
    assert anphon.Result("rel", "p").workdir.is_absolute()
    assert not any(hasattr(anphon.Result, n) for n in ("kappa", "scph", "qha"))


def test_prefix_follows_locate_tag(tmp_path):
    from anphon._api import _prefix_of

    def prefix(text):
        p = tmp_path / "x.in"
        p.write_text(text)
        return _prefix_of(p)

    good = "&general\n  MODE = phonons; prefix = si # comment\n/\n"
    assert prefix(good) == "si"
    # "#&general" is not the token "&general"; the real block comes later
    assert prefix("#&general\n PREFIX = no\n/\n" + good) == "si"
    assert prefix("&GENERAL\n # PREFIX = no\n PREFIX=yes\n/\n") == "yes"
    with pytest.raises(ValueError):
        prefix("&general\n MODE = phonons\n/\n&cell\n PREFIX = no\n/\n")
    with pytest.raises(ValueError):
        prefix("#&general\n PREFIX = no\n/\n")


def test_integer_fields_and_values():
    from anphon import format_input

    text = format_input(
        general=dict(
            TMIN=np.float32(0.5), EPS=1e-5, FLAG=np.bool_(True), N=np.int64(3)
        ),
        kpoint=dict(KPMODE=2.0, mesh=[np.int32(4), 4.0, 4]),
    )
    assert "TMIN = 0.5\n" in text and "EPS = 1e-05\n" in text
    assert "FLAG = 1\n" in text and "N = 3\n" in text
    assert "&kpoint\n  2\n  4 4 4\n/" in text
    path = [("G", [0, 0, 0], "X", [0.5, 0.5, 0], 51.0)]
    assert "G 0 0 0 X 0.5 0.5 0 51\n" in format_input(kpoint=dict(KPMODE=1, path=path))
    for bad in (
        dict(KPMODE=2.5, mesh=[4, 4, 4]),
        dict(KPMODE=2, mesh=[4, 4.5, 4]),
        dict(KPMODE=1, path=[("G", [0, 0, 0], "X", [0.5, 0.5, 0], 51.5)]),
        dict(KPMODE=True, mesh=[4, 4, 4]),
    ):
        with pytest.raises(ValueError):
            format_input(kpoint=bad)
    for bad in (float("nan"), np.inf, "a;b", "a=b"):
        with pytest.raises(ValueError):
            format_input(general=dict(X=bad))


def test_read_h5_decodes_nested_strings(tmp_path):
    import h5py
    from anphon import read_h5

    p = tmp_path / "s.h5"
    with h5py.File(p, "w") as f:
        f["s2d"] = np.array([[b"a", b"b"], [b"c", b"d"]])
        f["s2d"].attrs["names"] = np.array([b"x", b"y"])
        f.attrs["title"] = b"t"
    d = read_h5(p)
    assert d["s2d"] == [["a", "b"], ["c", "d"]]
    assert d["s2d@names"] == ["x", "y"] and d["@title"] == "t"
