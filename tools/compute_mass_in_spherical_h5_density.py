#!/usr/bin/env python

#
# ./compute_mass_in_spherical_h5_density.py <input.h5>
#
# where <input.h5> is an HDF5 file with a density on a spherical grid
#

import argparse
import sys

import h5py
import numpy as np

my_description = """
Computes the total mass and the radial extents of a density field given on a
3D spherical grid in an HDF5 file. The file must contain the datasets:
 - r, theta, phi: 1D coordinate arrays (theta: polar angle, 0..pi;
                  phi: azimuthal angle, 0..2*pi);
 - density:       3D array indexed as [theta, phi, r].
The angular coordinates are cell centers. The radial coordinates are either
grid nodes (default; the mass is integrated with the trapezoidal rule), or
cell centers (--r-centering cells). Mass in solar masses is only meaningful
if the file is in CGS units."""
my_usage = "./compute_mass_in_spherical_h5_density.py input.h5"

M_SUN_CGS = 1.98841e33  # [g]


def cell_edges(x, lo=None, hi=None):
    """Edges of the cells with centers x; outer edges are lo and hi if given,
    otherwise they are extrapolated."""
    mid = 0.5 * (x[1:] + x[:-1])
    lo = x[0] - 0.5 * (x[1] - x[0]) if lo is None else lo
    hi = x[-1] + 0.5 * (x[-1] - x[-2]) if hi is None else hi
    return np.concatenate([[lo], mid, [hi]])


def main():
    parser = argparse.ArgumentParser(
        description=my_description, usage=my_usage,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ifile", help="input HDF5 file")
    parser.add_argument("--r-centering", choices=["nodes", "cells"],
                        default="nodes",
                        help="radial coordinates are grid nodes or cell "
                             "centers [%(default)s]")
    parser.add_argument("--threshold", type=float, default=0.0,
                        help="density above this value counts as non-zero "
                             "for the extents of the matter [%(default)g]")
    args = parser.parse_args()

    with h5py.File(args.ifile, "r") as f:
        for name in ("r", "theta", "phi", "density"):
            if name not in f:
                sys.exit("ERROR: no dataset '%s' in %s" % (name, args.ifile))
        r, theta, phi = f["r"][:], f["theta"][:], f["phi"][:]
        rho = f["density"][:]
    if rho.shape != (theta.size, phi.size, r.size):
        sys.exit("ERROR: density has shape %s, expected [theta, phi, r] = %s"
                 % (rho.shape, (theta.size, phi.size, r.size)))

    # solid angle of the angular cells: (cos(theta1) - cos(theta2)) * dphi
    dmu = -np.diff(np.cos(cell_edges(theta, 0.0, np.pi)))
    dphi = 2.0 * np.pi / phi.size
    # density integrated over the solid angle, as a function of radius
    rho_r = np.einsum("tpr,t->r", rho, dmu) * dphi

    # a) total mass
    if args.r_centering == "nodes":
        y = rho_r * r**2
        mass = np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(r))
        rmin, rmax = r[0], r[-1]
    else:
        edges = cell_edges(r)
        mass = np.sum(rho_r * np.diff(edges**3) / 3.0)
        rmin, rmax = edges[0], edges[-1]

    # b) extents: of the grid, and of the matter (non-zero density)
    nonzero = np.nonzero((rho > args.threshold).any(axis=(0, 1)))[0]

    print("file:        %s" % args.ifile)
    print("grid:        Ntheta = %d, Nphi = %d, Nr = %d"
          % (theta.size, phi.size, r.size))
    print("total mass:  %.6e (%.6e Msun, if in CGS)"
          % (mass, mass / M_SUN_CGS))
    print("grid extent: rmin = %.6e, rmax = %.6e" % (rmin, rmax))
    if nonzero.size > 0:
        print("matter:      rmin = %.6e, rmax = %.6e  (density > %g)"
              % (r[nonzero[0]], r[nonzero[-1]], args.threshold))
    else:
        print("matter:      density <= %g everywhere" % args.threshold)
    print("density:     min = %.6e, max = %.6e" % (rho.min(), rho.max()))


if __name__ == "__main__":
    main()
