#!/usr/bin/env python
"""Check the warning about an imprecise input structure (SYMMETRIZE_STRUCTURE = 0).

alm suggest on the tracked Si input must not warn; the same input with a lattice
vector stretched by 1e-4 must warn and recommend SYMMETRIZE_STRUCTURE = 1, and
must stay silent when SYMMETRIZE_STRUCTURE = 1 is given.

Run:  python test_symmetry_deviation.py [path/to/alm]
"""

import os
import re
import subprocess
import sys

WARNING = "Consider SYMMETRIZE_STRUCTURE = 1"


def run_suggest(almbin, text, name):
    with open(name + ".in", "w") as f:
        f.write(text)
    ret = subprocess.run([almbin, name + ".in"], capture_output=True, text=True)
    if ret.returncode != 0:
        print("alm failed for %s:\n%s" % (name, ret.stdout[-2000:]))
        return None
    return ret.stdout


def runtest(almbin, project_root):
    with open("%s/example/Si/si_alm.in" % project_root) as f:
        text = f.read()
    text = re.sub(r"&optimize.*?/\n", "", text, flags=re.S)
    text = re.sub(r"MODE\s*=\s*\w+", "MODE = suggest", text)
    perturbed = text.replace("  1.0 0.0 0.0 # a1", "  1.0001 0.0 0.0 # a1")
    assert perturbed != text
    symmetrized = perturbed.replace("KD = Si", "KD = Si; SYMMETRIZE_STRUCTURE = 1")

    cases = [
        ("si_ref", text, False),
        ("si_perturbed", perturbed, True),
        ("si_symmetrized", symmetrized, False),
    ]
    failed = False
    for name, inp, expect in cases:
        out = run_suggest(almbin, inp, name)
        if out is None:
            return 1
        found = WARNING in out
        print(
            "%-15s warning %-7s (expected %s)"
            % (
                name,
                "printed" if found else "absent",
                "printed" if expect else "absent",
            )
        )
        failed |= found != expect
    return 1 if failed else 0


if __name__ == "__main__":
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if len(sys.argv) > 1:
        almbin = os.path.abspath(sys.argv[1])
    else:
        candidates = [
            "%s/_build/alm/alm" % project_root,
            "%s/build/alm/alm" % project_root,
        ]
        almbin = next((c for c in candidates if os.path.exists(c)), candidates[0])

    workdir = "%s/test/si_symdev" % project_root
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)

    info = runtest(almbin, project_root)
    print("Symmetry deviation warning --> %s" % ("pass" if info == 0 else "failed"))
    sys.exit(info)
