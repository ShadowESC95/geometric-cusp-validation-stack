# Cusp benchmark validation

The benchmark snapshot contains **480 pytest tests and 63 example groups**.
The 60th example group was added in commit
`ab067d159aa92c3abef94b4eaf8cb69eaaaa3a33` on `main`. The example runner defines the coverage of this snapshot.

## Reproduce the snapshot

From the repository root, use Python 3.10 or later and GNU Make. If you already
have the project environment active, start with `make test`.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
make test && make examples
```

The expected final lines are `480 passed` (with elapsed time) and
`PASS: all 63 example groups`. The user run on Linux reported both. The
clean-checkout verification used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0,
SymPy 1.14.0 and pytest 9.2.0. For a run record, save the full commit ID,
dependency versions, and command output. The numerical examples include
explicit tolerances, and the scrambled quadrature comparisons use fixed seeds.

Use the [example guide](../examples/README.md) for individual commands and
derivations. `make full` runs the broader research campaign with its own
dependencies and outputs.

## Results and scope

| Check | Result | Scope |
| --- | --- | --- |
| One-loop geometry and orientation | The adjacent-propagator exceptional divisor carries the logarithmic residue. The generic unweighted strict-cusp blow-up gives zero. The weighted strict-cusp chart recovers a nonzero residue and a variable cross-ratio. | Exact in the stated charts and orientations. |
| Finite CD weight and contour | Independent quadrature reproduces the finite CD function on the tested sheet and the contour mean $`2\pi^2/3`$. The capped weighted cycle has oriented quotient factor $`+1`$. | The benchmark normalization gives $`\Gamma^{(1)}=4`$ and $`\Gamma^{(2)}=-4\pi^2/3`$. With the published $`F^{(2)}`$ the same rule gives $`\Gamma^{(3)}=44\pi^4/45`$. |
| Lorentzian prescription | Chart transport, energy residues, finite-frame integrals, and regulator limits agree in the cases exercised by the suite. | The checks use the specified prescription, frame, and limit paths. |
| Regulated complex cap | The fixed auxiliary lift has regulated poles in the lower region. The diagonal cap with $`R=0.06`$ has a uniform denominator gap for $`x\ge2500\epsilon`$ when $`0<\epsilon\le10^{-4}`$. | The bound applies to the diagonal product cap. The transverse rim has additional degeneracies. |
| Integrated one-meridian band | At $`t=y/x=1`$, the lower band has a nonzero counterclockwise coefficient per unit $`t`$, approximately $`-11.2420175678-7.4161239685i`$ times $`\epsilon`$. | Example 60 integrates over the auxiliary triangle and the full band in the scaled $`x`$ direction on this slice. |

| Lorentzian edge chain | On the edge segment the Lorentzian contour in the $`\log z`$ plane deforms to the unit circle; off the segment it is contractible. | Local statement, exact identities plus numerics at $`10^{-7}`$. See [GATE_6_CLOSURE.md](GATE_6_CLOSURE.md) for what is argued and not proved. |

The benchmark suite is complete at this snapshot. A chain-level comparison
of the geometric and Lorentzian integration chains is given in
[GATE_6_CLOSURE.md](GATE_6_CLOSURE.md) and
[GATE_6_BOUNDS.md](GATE_6_BOUNDS.md), with the integration contour in
[GATE_6_CONTOUR.md](GATE_6_CONTOUR.md). It is unreviewed; the table in section 6
of the latter lists which statements are proved, computed, or argued. The
[Gate 6 proof document](GATE_6_PROOF.md) records the mathematical claims,
the transverse rim collision, and the two outstanding proof obligations.
