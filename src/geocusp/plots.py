from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


def save(fig, path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_amplituhedron_3d(points, hull, path, title='A_{1,n,4} affine slice'):
    fig = plt.figure(figsize=(8,7)); ax=fig.add_subplot(111, projection='3d')
    p=points
    ax.scatter(p[:,0],p[:,1],p[:,2],s=5,alpha=0.25)
    if hull is not None:
        for tri in hull.simplices:
            tri=tri[[0,1,2,0]]
            ax.plot(p[tri,0],p[tri,1],p[tri,2],linewidth=0.6,alpha=0.65)
    ax.set_xlabel('$Y_1/Y_0$'); ax.set_ylabel('$Y_2/Y_0$'); ax.set_zlabel('$Y_3/Y_0$')
    ax.set_title(title); save(fig,path)


def plot_3d_pca(P, path, title='G_+(k,n) pushforward in Plücker-PCA space'):
    X=P-P.mean(axis=0)
    U,S,Vt=np.linalg.svd(X, full_matrices=False)
    Q=X@Vt[:3].T
    fig=plt.figure(figsize=(8,7)); ax=fig.add_subplot(111,projection='3d')
    ax.scatter(Q[:,0],Q[:,1],Q[:,2],s=5,alpha=0.25)
    ax.set_xlabel('PC1'); ax.set_ylabel('PC2'); ax.set_zlabel('PC3'); ax.set_title(title)
    save(fig,path)


def plot_positivity(margins, path):
    fig,ax=plt.subplots(figsize=(8,5)); ax.hist(margins,bins=60); ax.set_xlabel('minimum ordered k-minor'); ax.set_ylabel('count'); ax.set_title('Positive-Grassmannian margin')
    save(fig,path)


def plot_complex_curve(curve, path, title='Complex image / linking scan'):
    fig,ax=plt.subplots(figsize=(7,7)); ax.plot(curve.real,curve.imag,linewidth=1.4); ax.axhline(0,linewidth=0.6); ax.axvline(0,linewidth=0.6); ax.set_aspect('equal','box'); ax.set_xlabel('Re'); ax.set_ylabel('Im'); ax.set_title(title); save(fig,path)


def plot_log_scaling(rhos, vals, path, reference_slope=-1.0, title='Radial scaling'):
    fig,ax=plt.subplots(figsize=(8,5)); ax.loglog(rhos,np.abs(vals),'o-',label='measured');
    if len(rhos)>=2:
        c=np.abs(vals[-1])/(rhos[-1]**reference_slope)
        ax.loglog(rhos,c*rhos**reference_slope,'--',label=f'ref. slope {reference_slope:g}')
    ax.set_xlabel(r'$\rho$'); ax.set_ylabel(r'$|f(\rho)|$'); ax.set_title(title); ax.legend(); save(fig,path)


def plot_qmc_convergence(ns, est, err, path):
    fig,ax=plt.subplots(figsize=(8,5)); ax.errorbar(ns,est,yerr=err,fmt='o-'); ax.set_xscale('log'); ax.set_xlabel('samples'); ax.set_ylabel('estimate'); ax.set_title('QMC / randomized digital-net convergence'); save(fig,path)


def plot_matrix_trajectory(step, comm, spec_path, path):
    fig,ax=plt.subplots(figsize=(8,5)); ax.plot(step,comm); ax.set_xlabel('measurement step'); ax.set_ylabel(r'$\sum_{\mu<\nu} -\mathrm{Tr}[X_\mu,X_\nu]^2/N$'); ax.set_title('Matrix-substrate Langevin trajectory'); save(fig,path)
    fig2,ax2=plt.subplots(figsize=(8,5)); ax2.hist(spec_path, bins=60); ax2.set_xlabel('eigenvalue'); ax2.set_ylabel('count'); ax2.set_title('Matrix eigenvalue distribution'); return fig2
