"""Silicon with the anphon Python API (ALAMODE tutorial).

Run from any directory; the outputs go to bands/, dos/, and kappa/ there:

    python /path/to/alamode/example/Si/python/si_python.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import anphon

REF = Path(__file__).resolve().parent.parent / "reference"  # example/Si/reference

# --- Step 1: phonon dispersion and DOS --------------------------------------
fc2 = str(REF / "si222.h5")  # absolute path: anphon runs inside workdir

bands = anphon.run(
    general=dict(PREFIX="si222", MODE="phonons", FCSFILE=fc2),
    kpoint=dict(
        KPMODE=1,
        path=[
            ("G", [0, 0, 0], "X", [0.5, 0.5, 0], 51),
            ("X", [0.5, 0.5, 1], "G", [0, 0, 0], 51),
            ("G", [0, 0, 0], "L", [0.5, 0.5, 0.5], 51),
        ],
    ),
    workdir="bands",
)
dos = anphon.run(
    general=dict(PREFIX="si222", MODE="phonons", FCSFILE=fc2, DELTA_E=2),
    kpoint=dict(KPMODE=2, mesh=[20, 20, 20]),
    workdir="dos",
)
print(bands.files(), dos.files())

# .bands: column 0 is the path coordinate, then 6 branches (cm^-1)
b = np.loadtxt(bands.path(".bands"))
d = np.loadtxt(dos.path(".dos"))  # column 0: frequency (cm^-1), column 1: total DOS
with open(bands.path(".bands")) as f:  # first two lines: labels and their positions
    labels = f.readline().split()[1:]
    ticks = [float(x) for x in f.readline().split()[1:]]

fig, (ax1, ax2) = plt.subplots(
    1,
    2,
    sharey=True,
    figsize=(6, 3.5),
    gridspec_kw=dict(width_ratios=[3, 1], wspace=0.05),
)
ax1.plot(b[:, 0], b[:, 1:], color="C0")
ax1.set_xticks(ticks, [r"$\Gamma$" if s == "G" else s for s in labels])
ax1.set_xlim(ticks[0], ticks[-1])
for t in ticks[1:-1]:
    ax1.axvline(t, color="gray", lw=0.5)
ax1.set_ylabel("Frequency (cm$^{-1}$)")
ax1.set_ylim(0, None)
ax2.plot(d[:, 1], d[:, 0], color="C0")
ax2.set_xlabel("DOS")
ax2.set_xticks([])
ax2.set_xlim(0, None)
fig.savefig("si_phonons.png", dpi=120, bbox_inches="tight")

# --- Step 2: frequencies and eigenvectors at arbitrary q --------------------
ph = anphon.Phonon(fc2)
q = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0.5, 0.5]]  # Gamma, X, L (fractional)
w, evec = ph.frequencies(q, eigenvectors=True)
print(w.shape, evec.shape)  # (3, 6) (3, 6, 6)
print(np.round(w, 1))  # cm^-1

# evec[iq, imode, 3 * iatom + ixyz]: the highest mode at Gamma (optical)
# moves the two atoms in opposite directions.
print(np.round(evec[0, 5].reshape(2, 3), 3))

# --- Step 3: thermal conductivity on a coarse mesh --------------------------
kdir = Path("kappa")
kdir.mkdir(exist_ok=True)
(kdir / "si222.kappa.h5").unlink(missing_ok=True)  # else anphon restarts from it

res = anphon.run(
    general=dict(
        PREFIX="si222",
        MODE="kappa",
        FCSFILE=str(REF / "si222_cubic.h5"),
        TMIN=100,
        TMAX=1000,
        DT=50,
    ),
    kpoint=dict(KPMODE=2, mesh=[8, 8, 8]),
    workdir=kdir,
)
k = res.h5("kappa")  # the raw contents of si222.kappa.h5
T = k["metadata"]["temperatures"]  # K
kap = k["kappa"]["kappa_peierls"]  # W/mK, shape (nT, 3, 3)
for Ti, ki in zip(T, kap):
    if Ti % 100 == 0:
        print(f"{Ti:6.0f} K {ki[0, 0]:8.1f} W/mK")

fig, ax = plt.subplots(figsize=(4, 3))
ax.plot(T, kap[:, 0, 0], "o-")
ax.set_xlabel("Temperature (K)")
ax.set_ylabel(r"$\kappa$ (W/mK)")
fig.savefig("si_kappa.png", dpi=120, bbox_inches="tight")
