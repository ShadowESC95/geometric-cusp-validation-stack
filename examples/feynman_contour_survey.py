"""An explicit integration contour for the frozen-loop integral.

See docs/GATE_6_CONTOUR.md. The real domain is displaced by x0 -> x0 + i V
with V built from four frame fields and two far fields. The example checks
the projection inequality behind the construction and surveys the contour on
a stress sample: correct side of every cone, small excursion of the phases,
and zeros only on the edges.
"""

import numpy as np

from geocusp import feynman_contour as fc


def check_feynman_contour():
    contour = fc.Contour()
    rng = np.random.default_rng(0)
    vertices = np.array([[0.0, 0.0], [0.0, -1.0], [1.0, 0.0]])
    for _ in range(5000):
        P = rng.normal(size=2) * 2
        d = P - fc.project_to_triangle(P)
        assert min(d @ (P - v) - d @ d for v in vertices) > -1e-12

    worst = -np.inf
    for _ in range(2000):
        x0 = contour.centre + rng.normal(size=4) * 3
        for j in range(4):
            v, bound = contour.frame_field(j, x0), contour.frame_bound(j, x0)
            for i in ((j - 1) % 4, j, (j + 1) % 4):
                worst = max(worst, fc.mdot(v, x0 - contour.x[i]) - bound)
    assert worst < 1e-9
    print(f"Lemma 5: V_j.(x0 - x_i) exceeds its bound by at most {max(worst, 0):.1e}")

    result = fc.survey(contour, 6000, seed=3)
    assert result["max_im_near"] <= 0 and result["max_excursion"] < 0.01
    assert result["min_ratio"] > 1e-9
    for region in ("--", "+-", "-+"):
        cert = fc.certificate(fc.Contour(region=region))
        assert cert["all"], [k for k, v in cert.items() if isinstance(v, bool) and not v]
        print(f"signs {region}: every condition of the proof holds; "
              f"||dV|| <= {cert['jacobian_bound']:.3f}, excursion <= {cert['excursion_bound']:.3f} rad")
    print(f"{result['near_count']} near-cone pairs: largest Im x0i^2 = {result['max_im_near']:.1e}")
    print(f"largest excursion of arg x0i^2 beyond [-pi, 0]: {result['max_excursion']:.1e} rad")
    print(f"smallest |x0i^2| / min(1, dist_to_edges^2): {result['min_ratio']:.1e}")
    bounds = [fc.edge_lower_bounds(contour, 2000, seed=k, cusp=k) for k in range(4)]
    assert min(bounds) > 0.5
    print("near-cusp lower bounds in units of nu:", ", ".join(f"{b:.0f}" for b in bounds))
    print("PASS: explicit contour stays on the Feynman side of every cone and meets them only on the edges")


if __name__ == "__main__":
    check_feynman_contour()
