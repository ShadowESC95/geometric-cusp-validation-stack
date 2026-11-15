# pySecDec handoff

Use pySecDec only after the geometric derivation has produced a regulated Feynman/parameter integral.
The intended workflow is:

1. derive the exact parameter representation in `src/geocusp`;
2. export a sector-decomposition input;
3. compute the Laurent expansion in the regulator;
4. identify the coefficient that corresponds to the cusp observable;
5. compare against the geometric period.

Do not use sector decomposition to define the missing geometric residue prescription.
