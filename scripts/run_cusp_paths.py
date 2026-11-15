import matplotlib.pyplot as plt
from geocusp.cusp_paths import strict_cusp_smoothing
from pathlib import Path
Path('figures').mkdir(exist_ok=True)
fig,ax=plt.subplots(figsize=(8,5))
for mode in ['symmetric','A_normal','B_normal']:
    rho,z=strict_cusp_smoothing(mode=mode)
    ax.semilogx(rho, z.real,label=mode)
ax.set_xlabel(r'$\rho$'); ax.set_ylabel(r'Re$(z)$'); ax.set_title('Strict-cusp smoothing: approach-path dependence of z'); ax.legend(); fig.savefig('figures/strict_cusp_approach_paths.png',dpi=180,bbox_inches='tight')
print('wrote figures/strict_cusp_approach_paths.png')
