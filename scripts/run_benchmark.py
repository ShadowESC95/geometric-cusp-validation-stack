#!/usr/bin/env python
"""Run the established four-point contour benchmark and archive raw samples."""
import json
from pathlib import Path
import numpy as np
from geocusp.contours import contour_average, benchmark_F0, benchmark_F1, gamma_coeffs
from geocusp.benchmark import cusp_coefficient_from_mean

out=Path('results'); fig=Path('figures'); out.mkdir(exist_ok=True); fig.mkdir(exist_ok=True)
computed = {}
reference = gamma_coeffs()
for order, (name, fn) in enumerate([('F0',benchmark_F0),('F1',benchmark_F1)], start=1):
    avg,z,vals=contour_average(fn)
    print(f'{name}: contour average={avg}')
    coefficient = cusp_coefficient_from_mean(order, avg)
    np.testing.assert_allclose(coefficient, reference[order], rtol=1e-6, atol=1e-12)
    computed[order] = float(coefficient.real)
    print(f'Gamma_cusp[{order}]: computed={coefficient.real:.12g}, reference={reference[order]:.12g}')
    phi=np.unwrap(np.angle(z))
    np.savetxt(out/f'{name}_contour.dat',np.c_[phi,vals.real,vals.imag], header='phi Re Im')
(out/'cusp_benchmark_coefficients.json').write_text(json.dumps(reference,indent=2))
(out/'cusp_benchmark_computed.json').write_text(json.dumps(computed, indent=2))
