# Architecture

The stack is organized as an executable dependency graph.

```text
                  exact geometry
                       |
        +--------------+----------------+
        |              |                |
        v              v                v
  Grassmannian      twistor          matrix substrate
  / amplituhedron   brackets         Langevin probe
        |              |                |
        v              v                v
  pushforward      cusp strata       emergent fields
        |              |
        +------+-------+
               v
        blow-up / residues
               |
               v
       linking / relative cycle
               |
               v
           z pushforward
               |
               v
          contour period
               |
               v
         Gamma_cusp benchmark
               |
               v
            sigma_IR
               |
       +-------+--------+
       |                |
       v                v
 continuum PDE      quantum boundary
       |                |
 FEniCSx/Firedrake  Kwant/QuTiP
 deal.II/FiPy/      |
 COMSOL             v
       |           entropy probes
       v
Meep/OpenFOAM
       |
       +---------------------+
       |                     |
       v                     v
 VASP vacuum/VEV       Geant4 particle layer
```

## Backend principle

Each external code attacks a different mathematical/physical layer. Python is the orchestrator and provenance layer. The native backends are deliberately not made to share the same discretization, so agreement between independent implementations is a meaningful validation signal.
