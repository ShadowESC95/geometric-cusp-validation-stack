"""Run the complete set of cusp validation examples from the repository root."""

from .cd_integral import check_cd_benchmark, check_continuation, check_contour_and_tail
from .chart_transport import check_chart_matching, check_edge_residue, check_transport
from .closed_pinch_curve import check_closed_pinch_curve
from .energy_contour import check_energy_contour
from .geometry import check_complex_ray, check_local_geometry, check_lorentzian_poles
from .lower_null_endpoint import check_lower_null_endpoint
from .null_endpoint import check_null_endpoint
from .null_insertion_crossing import check_null_insertion_crossing
from .null_pinch import check_null_pinch
from .pinch_intersection import check_pinch_intersection
from .spatial_pinch import check_spatial_pinch
from .spatial_tail import check_spatial_tail
from .tube_energy_sum import check_tube_energy_sum


def main():
    checks = (
        ("Local residues, winding and two-loop scaling", check_local_geometry),
        ("Finite CD integral and fixed-z scaling region", check_cd_benchmark),
        ("Complex z contour and uniform parameter-tail bound", check_contour_and_tail),
        ("Complex-ray poles and bypass jump", check_complex_ray),
        ("Lorentzian poles, partial residue and propagator signs", check_lorentzian_poles),
        ("Finite-cutoff Lorentzian continuation", check_continuation),
        ("Weighted twistor and regulator transport", check_transport),
        ("External pole at the chart edge", check_edge_residue),
        ("Oriented chart matching and integrated patch residue", check_chart_matching),
        ("Five-propagator energy integral and transverse pinch", check_energy_contour),
        ("Merging pinch roots and a spatial patch", check_spatial_pinch),
        ("Intersecting timelike pinch surfaces", check_pinch_intersection),
        ("Closed timelike pinch curve", check_closed_pinch_curve),
        ("Complete energy contour over the closed-pinch tube", check_tube_energy_sum),
        ("Large spatial momentum tail", check_spatial_tail),
        ("External null-edge cylinder", check_null_pinch),
        ("External null edge crossing the insertion pole", check_null_insertion_crossing),
        ("External null-edge endpoint", check_null_endpoint),
        ("Lower-branch null-edge endpoint", check_lower_null_endpoint),
    )
    for index, (name, check) in enumerate(checks, start=1):
        print(f"\n[{index}/{len(checks)}] {name}", flush=True)
        check()
    print(f"\nPASS: all {len(checks)} example groups")


if __name__ == "__main__":
    main()
