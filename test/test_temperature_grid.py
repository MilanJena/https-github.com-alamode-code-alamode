#!/usr/bin/env python
"""Python readers count TMIN/TMAX/DT points like System::get_num_temperature_points."""

import ast
import os

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
readers = ["tools/analyzer/anphonio.py", "tools/interpolate_alamode.py"]

# (TMIN, TMAX, DT, expected count); the first two are integral ranges whose
# quotient lands just under an integer in binary64; the last snaps a half
# integer away from zero, as std::round does.
cases = [
    (280, 280.4, 0.1, 5),
    (0, 1, 0.1, 11),
    (0, 1000, 10, 101),
    (0, 0.8, 0.3, 3),
    (0, 1, 0.3, 4),
    (0, 50000000.5, 1, 50000002),
]

for path in readers:
    tree = ast.parse(open(os.path.join(root, path)).read())
    func = next(
        n
        for n in tree.body
        if isinstance(n, ast.FunctionDef) and n.name == "num_temperature_points"
    )
    namespace = {}
    exec(compile("import math\n" + ast.unparse(func), path, "exec"), namespace)
    for t_min, t_max, dt, expected in cases:
        got = namespace["num_temperature_points"](t_min, t_max, dt)
        assert got == expected, (
            f"{path}: ({t_min}, {t_max}, {dt}) gives {got}, expected {expected}"
        )
    print(f"{path} --> pass")
