"""Python layer over ``run_input_file``: input files from dicts, HDF5 result readers,
and harmonic phonons at arbitrary q (API stability: docs/source/anphon_python.rst).

Everything here is file based: ``run`` writes an ordinary anphon input file and the
results are read back from the HDF5 files anphon writes (docs/source/hdf5_format.rst).
"""

import copy
import glob
import math
import numbers
import os
import re
import tempfile
from pathlib import Path

import numpy as np

from ._anphon import abort_all, comm_size, run_input_file

__all__ = ["run", "Result", "Phonon", "format_input", "read_h5"]


# ---------------------------------------------------------------------------
# input text
# ---------------------------------------------------------------------------


def _items(v):
    """Flatten nested sequences / arrays into scalars (strings stay whole)."""
    if (
        isinstance(v, (str, bytes))
        or not hasattr(v, "__iter__")
        or getattr(v, "shape", None) == ()
    ):
        yield v
    else:
        for x in v:
            yield from _items(x)


def _scalar(x):
    if isinstance(x, bytes):
        x = x.decode()
    if isinstance(x, (bool, np.bool_)):
        return "1" if x else "0"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, (float, np.floating)):
        if not math.isfinite(x):
            raise ValueError(f"anphon: non-finite value {x!r}")
        return repr(float(x))
    return str(x)


def _int(x):
    """An integer-only field (KPMODE, mesh, npts): integral numbers only."""
    if isinstance(x, (bool, np.bool_)):
        pass
    elif isinstance(x, numbers.Integral):
        return int(x)
    elif isinstance(x, numbers.Real) and float(x).is_integer():
        return int(x)
    raise ValueError(f"anphon: expected an integer, got {x!r}")


def _fmt(v):
    """One input value: bools -> 1/0, sequences joined by spaces."""
    s = " ".join(_scalar(x) for x in _items(v))
    if any(c in s for c in "\n;#="):
        raise ValueError(f"anphon: value {s!r} contains a newline, ';', '#' or '='")
    return s


def _rows(rows):
    return [" ".join(_scalar(x) for x in _items(r)) for r in rows]


def _keyvals(d):
    return [f"{k.upper()} = {_fmt(v)}" for k, v in d.items() if v is not None]


def _cell(d):
    """&cell: dict(a=<factor>, lattice=<3x3, rows a1 a2 a3>) or a bare 3x3 lattice."""
    if not isinstance(d, dict):
        d = {"lattice": d}
    lat = np.asarray(d["lattice"], dtype=float)
    if lat.shape != (3, 3):
        raise ValueError("anphon: cell lattice must be 3x3 (rows a1, a2, a3)")
    return [_scalar(d.get("a", 1.0))] + _rows(lat)


def _kpoint(d):
    """&kpoint: KPMODE line, then mesh (2) / path (1) / qpoints (0) / raw lines."""
    mode = _int(d["KPMODE"] if "KPMODE" in d else d["kpmode"])
    if "lines" in d:
        body = [str(x) for x in d["lines"]]
    elif mode == 2:
        body = _rows([[_int(n) for n in d["mesh"]]])
    elif mode == 1:  # (label1, q1, label2, q2, npts) per segment
        body = _rows([[l1, q1, l2, q2, _int(n)] for l1, q1, l2, q2, n in d["path"]])
    else:
        body = _rows(np.reshape(np.asarray(d["qpoints"], dtype=float), (-1, 3)))
    return [str(mode)] + body


# positional sections; every other section is "KEY = value" lines
_FORMATTERS = {"cell": _cell, "kpoint": _kpoint}


def format_input(**sections):
    """Return anphon input text for the given sections (``general=dict(...)``, ...).

    A section value is
    - a ``str``: written verbatim as the body of the block (escape hatch);
    - ``cell``: ``dict(a=..., lattice=[[a1], [a2], [a3]])`` or a bare 3x3 lattice;
    - ``kpoint``: ``dict(KPMODE=2, mesh=[n1, n2, n3])``,
      ``dict(KPMODE=1, path=[("G", [0, 0, 0], "X", [.5, .5, 0], 51), ...])``,
      ``dict(KPMODE=0, qpoints=[[qx, qy, qz], ...])`` or ``dict(KPMODE=k, lines=[...])``;
    - a ``dict`` (other sections): ``KEY = value`` lines; lists are joined by spaces,
      bools become 1/0, ``None`` entries are dropped;
    - a sequence of rows (e.g. ``strain`` as a 3x3 matrix, ``displace``): one line per row.
    Only well-formedness is checked here; anphon's parser validates the content.
    """
    out = []
    for name, body in sections.items():
        if body is None:
            continue
        if isinstance(body, str):
            lines = body.strip("\n").splitlines()
        elif name in _FORMATTERS:
            try:
                lines = _FORMATTERS[name](body)
            except KeyError as exc:
                raise ValueError(f"anphon: section {name} lacks {exc}") from None
        elif isinstance(body, dict):
            lines = _keyvals(body)
        else:
            lines = _rows(body)
        out += [f"&{name}"] + [f"  {x}" for x in lines] + ["/", ""]
    return "\n".join(out)


def _write_atomic(path, text):
    # Under MPI every rank writes the same text: a private temp file + os.replace
    # makes each rank see a complete file without a barrier (also on NFS-like
    # shared filesystems), so no mpi4py dependency.
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            f.write(text)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def _prefix_of(input_path):
    """PREFIX as anphon reads it: InputParser::locate_tag finds the first
    whitespace-delimited token equal to "&general" (case-insensitive, comments are
    not stripped there), skips one character, then get_var_dict reads lines up to
    one starting with "/", dropping "#" comments and splitting on ";" and "="."""
    text = Path(input_path).read_text()
    tag = next(
        (m for m in re.finditer(r"\S+", text) if m.group().lower() == "&general"), None
    )
    if tag is not None:
        for line in text[tag.end() + 1 :].splitlines():
            line = line.split("#")[0].lstrip()
            if line.startswith("/"):
                break
            for entry in line.split(";"):
                k, _, v = entry.partition("=")
                if k.strip().upper() == "PREFIX":
                    return v.strip()
    raise ValueError(f"anphon: no PREFIX found in &general of {input_path}")


def run(input=None, *, workdir=None, **sections):
    """Run anphon and return a :class:`Result`.

    ``run("si.in", workdir=None)`` runs an existing input file (as ``anphon si.in``).
    ``run(general=dict(PREFIX="si", MODE="phonons", ...), kpoint=dict(...), ...,
    workdir="run1")`` writes ``workdir/<PREFIX>.in`` from the sections (see
    :func:`format_input`) and runs it. ``workdir`` (default: the cwd) is created if
    missing; anphon runs inside it, so relative file names in the input (FCSFILE, ...)
    and all outputs are relative to ``workdir``.

    The keyword-only argument names of ``run`` (``workdir``, and future ones such as
    ``comm``) are reserved; every other keyword is an input section.

    Errors: in a serial run, errors while preparing the input raise ordinary Python
    exceptions (TypeError for an input file plus sections, ValueError for a missing
    ``general``/PREFIX or a bad value, OSError from the filesystem), and errors that
    anphon itself reports raise RuntimeError. Under MPI with more than one process,
    both kinds abort all ranks (also a preparation error on one rank only).

    Collective under MPI (every rank must call it with the same arguments). Outputs
    are written by rank 0; every rank returns only after they are complete.
    """
    try:  # preparation: rank-local, may fail on one rank only
        wd, prefix, path = _prepare(input, workdir, sections)
    except Exception as exc:
        if comm_size() > 1:  # the other ranks would wait in the collective run
            abort_all(f"anphon.run: {exc}")
        raise
    run_input_file(path, wd)
    return Result(wd, prefix)


def _prepare(input, workdir, sections):
    wd = Path(workdir if workdir is not None else os.getcwd()).absolute()
    wd.mkdir(parents=True, exist_ok=True)
    if input is not None:
        if sections:
            raise TypeError(
                "anphon.run: give either an input file or sections, not both"
            )
        return wd, _prefix_of(input), Path(input).absolute()
    if not isinstance(sections.get("general"), dict) or "PREFIX" not in {
        k.upper() for k in sections["general"]
    }:
        raise ValueError(
            "anphon.run: general=dict(PREFIX=..., MODE=..., ...) is required"
        )
    prefix = next(
        _fmt(v) for k, v in sections["general"].items() if k.upper() == "PREFIX"
    )
    path = wd / f"{prefix}.in"
    _write_atomic(path, format_input(**sections))
    return wd, prefix, path


# ---------------------------------------------------------------------------
# results
# ---------------------------------------------------------------------------


def _attr(v):
    """bytes -> str, also inside (nested) string arrays, which become lists."""
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, np.ndarray) and v.dtype.kind in "SO":
        return _attr(v.tolist())
    if isinstance(v, list):
        return [_attr(x) for x in v]
    return v


def read_h5(path):
    """Read a whole HDF5 file into a nested dict.

    Groups become dicts and datasets numpy arrays (scalars as numpy scalars, strings
    as ``str``). Attributes are stored next to them under ``"@name"`` for an attribute
    of the group itself and ``"dataset@name"`` for an attribute of a dataset.
    """
    import h5py

    def walk(g):
        d = {f"@{k}": _attr(v) for k, v in g.attrs.items()}
        for name, obj in g.items():
            if isinstance(obj, h5py.Group):
                d[name] = walk(obj)
            else:
                val = obj[()]
                if isinstance(val, bytes) or (
                    isinstance(val, np.ndarray) and val.dtype.kind in "SO"
                ):
                    val = _attr(val)
                d[name] = val
                d.update({f"{name}@{k}": _attr(v) for k, v in obj.attrs.items()})
        return d

    with h5py.File(path, "r") as f:
        return walk(f)


class Result:
    """Outputs of one anphon run: ``workdir`` (absolute Path) and ``prefix``."""

    def __init__(self, workdir, prefix):
        self.workdir = Path(workdir).absolute()
        self.prefix = prefix

    def __repr__(self):
        return f"Result(workdir={str(self.workdir)!r}, prefix={self.prefix!r})"

    def path(self, suffix):
        """``workdir / (prefix + suffix)``, e.g. ``r.path(".bands")``."""
        return self.workdir / f"{self.prefix}{suffix}"

    def files(self):
        """Sorted names of the files ``<prefix>.*`` in ``workdir``."""
        return sorted(
            p.name for p in self.workdir.glob(f"{glob.escape(self.prefix)}.*")
        )

    def h5(self, kind):
        """Raw contents of ``<workdir>/<prefix>.<kind>.h5`` as read by :func:`read_h5`,
        e.g. ``r.h5("kappa")`` (also "scph", "qha", "selfenergy", "scph_fe", "evec", ...).

        The layout is the raw file schema, versioned separately in
        docs/source/hdf5_format.rst. Raises FileNotFoundError if the file is absent.
        (The names kappa(), scph(), qha() are kept free for typed results later.)
        """
        p = self.path(f".{kind}.h5")
        if not p.is_file():
            raise FileNotFoundError(f"anphon: {p} not found")
        return read_h5(p)


# ---------------------------------------------------------------------------
# harmonic phonons at arbitrary q
# ---------------------------------------------------------------------------

# &general tags that Phonon sets itself (or takes as constructor arguments)
_PHONON_OWN = {
    "PREFIX",
    "MODE",
    "FCSFILE",
    "BORNINFO",
    "NONANALYTIC",
    "FILE_FORMAT",
    "KD",
}


class Phonon:
    """Harmonic phonons at arbitrary q from a force-constant file (v1: file based).

    ``Phonon(fcsfile, *, cell=None, kd=None, borninfo=None, nonanalytic=0,
    general=None)``. Settings are validated and files loaded at construction or at
    the first calculation (an implementation detail that may change). In v1 each
    :meth:`frequencies` call runs ``MODE = phonons`` with ``KPMODE = 0`` in a
    temporary directory and reads ``<prefix>.evec.h5``.

    ``cell`` (``dict(a=..., lattice=...)`` or a 3x3 lattice in bohr, as ``&cell``) is
    needed only when the FCSFILE carries no primitive cell (an alm XML file does not,
    an HDF5 file from alm 2.0 does); ``kd`` (list of element names) is optional likewise.
    ``general`` is a dict of extra ``&general`` tags (MASS, TOLERANCE, VERBOSITY, ...);
    VERBOSITY defaults to 0 here, which prints nothing but warnings;
    the tags Phonon sets itself (PREFIX, MODE, FCSFILE, BORNINFO, NONANALYTIC,
    FILE_FORMAT, KD) raise ValueError there, in any letter case. The settings are
    copied, so later changes to the caller's objects have no effect.

    For ``nonanalytic > 0``, KPMODE = 0 takes the direction of the non-analytic term as
    zero at reciprocal-lattice Gamma points, i.e. no LO-TO splitting exactly at Gamma.

    Serial only in v1: under MPI each rank would get a different temporary directory
    and only rank 0 writes the output, so it raises RuntimeError when MPI_COMM_WORLD
    has more than one process.
    """

    def __init__(
        self,
        fcsfile,
        *,
        cell=None,
        kd=None,
        borninfo=None,
        nonanalytic=0,
        general=None,
    ):
        if general is None:
            general = {}
        elif not isinstance(general, dict):
            raise TypeError("anphon.Phonon: general must be a dict of &general tags")
        general = copy.deepcopy(general)
        own = sorted(k for k in general if str(k).upper() in _PHONON_OWN)
        if own:
            raise ValueError(f"anphon.Phonon sets {own} itself; use the arguments")
        self._fcsfile = os.path.abspath(fcsfile)
        self._cell = copy.deepcopy(cell)
        self._kd = copy.deepcopy(kd)
        self._borninfo = os.path.abspath(borninfo) if borninfo else None
        self._nonanalytic = nonanalytic
        self._general = general

    @property
    def fcsfile(self):
        """Absolute path of the force-constant file."""
        return self._fcsfile

    def frequencies(self, qpoints_frac, eigenvectors=False):
        """Frequencies (cm^-1, shape (nq, nmodes)) at q points in fractional units of
        the reciprocal primitive lattice, and with ``eigenvectors=True`` also the
        complex eigenvectors of the dynamical matrix, shape (nq, nmodes, nmodes):
        ``evec[iq, imode, 3 * iatom + ixyz]`` (axis 1 is the mode, as in
        ``/Eigenvalues/polarization_vectors``). Imaginary modes come out negative.
        ``qpoints_frac`` must have shape (3,) or (nq, 3); anything else raises ValueError.
        """
        if comm_size() > 1:
            raise RuntimeError("anphon.Phonon is serial only in v1")
        q = np.asarray(qpoints_frac, dtype=float)
        if q.shape == (3,):
            q = q[None, :]
        elif q.ndim != 2 or q.shape[1] != 3:
            raise ValueError(
                f"anphon.Phonon: q points must have shape (3,) or (nq, 3), got {q.shape}"
            )
        general = dict(
            PREFIX="phonon",
            MODE="phonons",
            FCSFILE=self._fcsfile,
            KD=self._kd,
            BORNINFO=self._borninfo,
            NONANALYTIC=self._nonanalytic,
            FILE_FORMAT="h5",
        )
        # quiet by default; anphon.run keeps the CLI default
        if not any(str(k).upper() == "VERBOSITY" for k in self._general):
            general["VERBOSITY"] = 0
        general.update(self._general)
        import h5py

        with tempfile.TemporaryDirectory(prefix="anphon_phonon_") as tmp:
            r = run(
                general=general,
                cell=self._cell,
                analysis=dict(PRINTEVEC=1),
                kpoint=dict(KPMODE=0, qpoints=q),
                workdir=tmp,
            )
            with h5py.File(r.path(".evec.h5"), "r") as f:
                freq = np.array(f["Eigenvalues/frequencies"])
                if eigenvectors:
                    p = f["Eigenvalues/polarization_vectors"][()]
                    evec = p[..., 0] + 1j * p[..., 1]
        return (freq, evec) if eigenvectors else freq
