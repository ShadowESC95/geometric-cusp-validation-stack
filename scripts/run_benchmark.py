#!/usr/bin/env python
"""Run the established four-point contour benchmark and archive raw samples."""
import json
from pathlib import Path
import numpy as np
from geocusp.contours import contour_average, benchmark_F0, benchmark_F1, gamma_coeffs

out=Path('results'); fig=Path('figures'); out.mkdir(exist_ok=True); fig.mkdir(exist_ok=True)
for name,fn in [('F0',benchmark_F0),('F1',benchmark_F1)]:
    avg,z,vals=contour_average(fn)
    print(f'{name}: contour average={avg}')
    phi=np.unwrap(np.angle(z))
    np.savetxt(out/f'{name}_contour.dat',np.c_[phi,vals.real,vals.imag], header='phi Re Im')
(out/'cusp_benchmark_coefficients.json').write_text(json.dumps(gamma_coeffs(),indent=2))
