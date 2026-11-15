from geocusp.linking import exceptional_ratio_curve,winding,pullback_dlog_check
from geocusp.visualization import residue_map_plot,linking_plot
eta,z,phi=exceptional_ratio_curve(); residue_map_plot(eta,z,'figures/exceptionsal_divisor_z_map.png'); linking_plot(z,'figures/linking_cycle_z_plane.png',winding(z))
print('winding=',winding(z)); print('pullback_dlog_error,winding=',pullback_dlog_check())
