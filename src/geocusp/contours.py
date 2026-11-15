from __future__ import annotations
import numpy as np
import mpmath as mp


def contour_points(radius=1.0, n=4096, phase0=-np.pi):
    phi = np.linspace(phase0, phase0+2*np.pi, n, endpoint=False)
    return radius*np.exp(1j*phi), phi


def contour_average(f, radius=1.0, n=4096):
    z, phi = contour_points(radius, n)
    vals = np.asarray([f(complex(zz)) for zz in z], dtype=complex)
    # (1/(2*pi*i)) ∮ dz/z f(z) = mean_theta f(exp(i theta)).
    avg = np.mean(vals)
    return avg, z, vals


def winding_number(curve: np.ndarray, point: complex=0j) -> float:
    shifted = np.asarray(curve) - point
    closed = np.append(shifted, shifted[0])
    arg = np.unwrap(np.angle(closed))
    return float((arg[-1] - arg[0]) / (2 * np.pi))


def benchmark_F0(z):
    return -1.0 + 0.0j


def benchmark_F1(z):
    return np.log(complex(z))**2 + np.pi**2


def gamma_weak(g):
    return 4*g**2 - (4*np.pi**2/3)*g**4 + (44*np.pi**4/45)*g**6


def gamma_coeffs():
    return {1:4.0, 2:-4*np.pi**2/3, 3:44*np.pi**4/45}
