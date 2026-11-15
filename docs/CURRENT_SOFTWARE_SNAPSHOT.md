# Current software snapshot

Checked against project documentation on 2026-09-18. The stack avoids pinning native/HPC packages because site toolchains differ; install the version supported by the user's machine and record it in `results/manifest.json`.

- FEniCSx/DOLFINx: official documentation currently exposes DOLFINx 0.11.0.post0 and UFL 2026.1.0, with development `main` documentation also available.
- Firedrake: current documentation exposes an actively tested development `main` branch; stable installation should follow the site's pinned environment.
- deal.II: 9.8.0 is current in the project news feed as of 2026-08-09.
- FiPy: current documentation is 4.0.3+development snapshot and describes a finite-volume PDE framework.
- QuTiP: current documentation lists 5.3.1 as the latest release.
- OpenFOAM: version 14 was released 2026-07-14, with patch release 20260724 documented in July 2026.
- pySecDec: current project repository describes sector-decomposition evaluation of dimensionally regulated parameter integrals.

The repository's backend adapters are version-aware in documentation but deliberately not hard-pinned to a single native installation.
