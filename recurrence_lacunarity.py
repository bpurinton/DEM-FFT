"""
Normalized recurrence lacunarity (nRL) of spatial recurrence plots.

Implements the measure described in:

* Bielig, C., Braun, T., Palmero, M. S., Rheinwalt, A., and Marwan, N.:
  Recurrence lacunarity of spatial data, Eur. Phys. J. Spec. Top., 2026.
  https://doi.org/10.1140/epjs/s11734-026-02569-4

* Braun, T., Unni, V. R., Sujith, R. I., Kurths, J., and Marwan, N.:
  Detection of dynamical regime transitions with lacunarity as a multiscale
  recurrence quantification measure, Nonlinear Dynamics, 2021.
  https://doi.org/10.1007/s11071-021-06457-5

For an ny x nx image the spatial recurrence plot (RP) is a 4D binary array
R[y1, x1, y2, x2] = 1 if |z[y1, x1] - z[y2, x2]| <= eps. It has (ny*nx)**2
entries, so it is never stored here: it is built one band of image rows at a
time and immediately reduced to box masses.
"""
import numpy as np


def recurrence_threshold(z, rr=0.2, tol=1e-4):
    """
    Find the threshold eps that gives a recurrence rate of `rr`, i.e. the
    fraction of all pixel pairs (including self-pairs) with |zi - zj| <= eps.

    For scalar data the pair count for a given eps is exact from the sorted
    values, so we simply bisect on eps instead of sampling pair distances.
    """
    s = np.sort(z[np.isfinite(z)].ravel())
    n = s.size

    def rate(eps):
        hi = np.searchsorted(s, s + eps, side='right')
        lo = np.searchsorted(s, s - eps, side='left')
        return np.sum(hi - lo) / n**2

    lo, hi = 0.0, s[-1] - s[0]
    while hi - lo > tol * (s[-1] - s[0]):
        mid = (lo + hi) / 2
        if rate(mid) < rr:
            lo = mid
        else:
            hi = mid
    return hi, rate(hi)


def box_masses(z, eps, w, max_chunk_bytes=2**28):
    """
    Masses of the non-overlapping 4D boxes of size w x w x w x w in the spatial
    RP of `z`, anchored at the upper-left corner (0, 0, 0, 0). Pixels that do
    not fill a whole box along the right/bottom edge are dropped.

    Returns an array of shape (Py, Px, Py, Px): the number of recurrences
    between every pair of w x w image patches.
    """
    Py, Px = z.shape[0] // w, z.shape[1] // w
    z = z[:Py * w, :Px * w].astype(np.float32)
    zf = z.ravel()

    # process one row of patches (w image rows) at a time, splitting further
    # into sub-bands of patch columns so the float32 distances fit in memory
    cols_per_chunk = max(1, int(max_chunk_bytes // (4 * w * w * zf.size)))
    M = np.empty((Py, Px, Py, Px), dtype=np.int32)
    for py in range(Py):
        band = z[py * w:(py + 1) * w]  # (w, Px*w)
        for px0 in range(0, Px, cols_per_chunk):
            px1 = min(Px, px0 + cols_per_chunk)
            vals = band[:, px0 * w:px1 * w]  # (w, n*w)
            d = np.subtract(vals[:, :, None], zf[None, None, :])
            R = np.abs(d, out=d) <= eps
            # sum each (w x w) source patch against each (w x w) target patch
            R = R.reshape(w, px1 - px0, w, Py, w, Px, w)
            M[py, px0:px1] = R.sum(axis=(0, 2, 4, 6))
    return M


def nrl(z, box_sizes, rr=0.2, eps=None, verbose=False):
    """
    Normalized recurrence lacunarity curve of the spatial RP of `z`.

    Lambda(w) = <M^2> / <M>^2 over all boxes of size w (Eq. 4), computed for
    both the RP and its inverse (masses w**4 - M), and combined as

        nRL(w) = 2 - (1 / Lambda(w) + 1 / Lambda_inv(w))      (Eq. 5)

    which is ~0 for a homogeneous RP and ~1 for a heterogeneous one.
    """
    if eps is None:
        eps, rr = recurrence_threshold(z, rr)

    # building the RP is the expensive part, so build the w = 2 masses once
    # and get every even w by summing k x k x k x k blocks of them
    M2 = box_masses(z, eps, 2) if any(w % 2 == 0 for w in box_sizes) else None
    out = []
    for w in box_sizes:
        if w % 2 == 0:
            k = w // 2
            Py, Px = z.shape[0] // w, z.shape[1] // w
            M = M2[:Py * k, :Px * k, :Py * k, :Px * k]
            M = M.reshape(Py, k, Px, k, Py, k, Px, k).sum(axis=(1, 3, 5, 7))
        else:
            M = box_masses(z, eps, w)
        M = M.astype(np.float64).ravel()
        Minv = w**4 - M
        lam = np.mean(M**2) / np.mean(M)**2
        lam_inv = np.mean(Minv**2) / np.mean(Minv)**2
        out.append(2 - (1 / lam + 1 / lam_inv))
        if verbose:
            print("w = {:3d}: nRL = {:.4f}".format(w, out[-1]))
    return np.array(out), eps, rr


def radial_spectrum(z, step):
    """
    Radially collapsed 1D DFT periodogram of `z`, following the same steps as
    example_analysis.ipynb (plane detrend, Hann window, power-of-2 padding).
    Returns frequency (1/m) and power, sorted by frequency.
    """
    ny, nx = z.shape
    x, y = np.meshgrid(range(nx), range(ny))
    A = np.vstack([x.ravel(), y.ravel(), np.ones(x.size)]).T
    fit = np.linalg.lstsq(A, z.ravel(), rcond=None)[0]
    z = z - (fit[0] * x + fit[1] * y + fit[2])

    hann_2d = np.sqrt(np.outer(np.hanning(ny), np.hanning(nx)))
    hann_weight = np.sum(hann_2d ** 2)
    z = z * hann_2d

    L = int(2**(np.ceil(np.log(max(nx, ny)) / np.log(2))))
    fft_shift = np.fft.fftshift(np.fft.fftn(z, (L, L)))
    c = L // 2
    fft_shift[c, c] = 0
    p2d = np.abs(fft_shift)**2 / (L * L * hann_weight)

    kx, ky = np.meshgrid(np.arange(L) - c, np.arange(L) - c)
    f2d = np.sqrt(kx**2 + ky**2) / (L * step)

    # keep the non-redundant half of the spectrum
    p1d = p2d[:, :c + 1].copy()
    f1d = f2d[:, :c + 1].copy()
    f1d[c:, c] = -1
    f1d, p1d = f1d.ravel(), p1d.ravel()
    keep = f1d > 0
    I = np.argsort(f1d[keep])
    return f1d[keep][I], 2 * p1d[keep][I]
