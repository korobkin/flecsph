#!/usr/bin/env python
"""
Inspect a density field on a 3D spherical grid stored in an HDF5 file.

Usage:
  ./inspect_spherical_h5_density.py input.h5
  ./inspect_spherical_h5_density.py input.h5 --plot-profile
  ./inspect_spherical_h5_density.py input.h5 --output-radial-profile
  ./inspect_spherical_h5_density.py input.h5 --output-normalized-radial-profile

What it reports:
  - total mass;
  - radial extent of the grid (rmin, rmax);
  - radial extent of the matter, i.e. of the region with non-zero density;
  - minimum and maximum density.

With --plot-profile, it also:
  - verifies that the density is spherically symmetric: the radial profiles
    density[j, k, :] must be identical for all theta and phi;
  - plots the radial profile into <input>_profile.png in the current
    directory.

With --output-radial-profile, it also:
  - verifies the spherical symmetry, as above;
  - writes the radial profile into the text file <input>.dat in the current
    directory, with two columns: 1:radius[cm] 2:density[g/cm3].
    If the density is not spherically symmetric, the density column is the
    average over the solid angle.

With --output-normalized-radial-profile, it also:
  - verifies the spherical symmetry, as above;
  - writes the radial profile, rescaled to unit radius and unit total mass,
    into the text file <input>_normalized.dat in the current directory.
    This is the 4-column format of the FleCSPH 1D density profiles
    (density_profile = "from_file"):
      1: r        radius divided by R0, from 0 to 1; R0 is the last radius
                  of the grid
      2: rho      density divided by rho0 = M/R0^3, where M is the total mass
      3: mass     4*pi * integral of rho r^2 dr from 0 to r (in the rescaled
                  variables), which grows from 0 to 1
      4: drho/dr  derivative of column 2 with respect to column 1, computed
                  with centered finite differences
    The profile is treated as point values at the radii of the grid, and
    both M and column 3 are integrated with the trapezoidal rule. If the grid
    does not start at r = 0, the profile is continued to r = 0 with zero
    density.

Input file, required datasets:
  r        1D, radius
  theta    1D, polar angle, 0..pi
  phi      1D, azimuthal angle, 0..2*pi
  density  3D, indexed as [theta, phi, r]

Grid conventions:
  - theta and phi are cell centers;
  - r is grid nodes by default, and the mass is integrated in radius with
    the trapezoidal rule; use --r-centering cells if r is cell centers.

Units:
  All quantities are printed in the units of the file. The mass in solar
  masses is only meaningful if the file is in CGS units.
"""

import argparse
import os
import sys

import h5py
import numpy as np

M_SUN_CGS = 1.98841e33  # [g]


def cell_edges(x, lo=None, hi=None):
    """Edges of the cells with centers x; outer edges are lo and hi if given,
    otherwise they are extrapolated."""
    mid = 0.5 * (x[1:] + x[:-1])
    lo = x[0] - 0.5 * (x[1] - x[0]) if lo is None else lo
    hi = x[-1] + 0.5 * (x[-1] - x[-2]) if hi is None else hi
    return np.concatenate([[lo], mid, [hi]])


def check_symmetry(rho, r, theta, phi):
    """Verify that the radial profiles rho[j, k, :] are identical for all
    theta and phi. Returns True if they are."""
    spread = rho.max(axis=(0, 1)) - rho.min(axis=(0, 1))  # at each radius
    if spread.max() == 0.0:
        print("symmetry:    OK, radial profiles are identical for all "
              "%d x %d angles" % (theta.size, phi.size))
        return True
    i = spread.argmax()
    j, k = np.unravel_index(np.abs(rho[:, :, i] - rho[0, 0, i]).argmax(),
                            rho.shape[:2])
    print("symmetry:    WARNING: radial profiles are NOT identical")
    print("             max difference = %.6e (%.3e of the max density)"
          % (spread.max(), spread.max() / np.abs(rho).max()))
    print("             at r = %.6e, between (theta, phi) = (%g, %g) "
          "and (%g, %g)" % (r[i], theta[0], phi[0], theta[j], phi[k]))
    return False


def plot_profile(ifile, rho, r, symmetric):
    """Plot the radial density profile into <ifile>_profile.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    if symmetric:
        ax.plot(r, rho[0, 0, :], lw=1.8, label="density")
    else:
        # not symmetric: show the spread over the angles
        ax.fill_between(r, rho.min(axis=(0, 1)), rho.max(axis=(0, 1)),
                        alpha=0.3, label="min .. max over angles")
        ax.plot(r, rho.mean(axis=(0, 1)), lw=1.8, label="mean over angles")
        ax.legend()
    ax.set_xlabel("r")
    ax.set_ylabel("density")
    ax.set_title(os.path.basename(ifile), fontsize=9)
    ax.grid(color="0.9")
    fig.tight_layout()
    ofile = os.path.splitext(os.path.basename(ifile))[0] + "_profile.png"
    fig.savefig(ofile, dpi=150)
    print("profile:     written to %s" % ofile)


def write_profile(ifile, r, profile):
    """Write the radial density profile into the text file <ifile>.dat."""
    ofile = os.path.splitext(os.path.basename(ifile))[0] + ".dat"
    np.savetxt(ofile, np.column_stack((r, profile)), fmt="%.12e",
               header="1:radius[cm] 2:density[g/cm3]")
    print("profile:     written to %s" % ofile)


def write_normalized_profile(ifile, r, profile):
    """Write the radial profile rescaled to unit radius and unit total mass
    into the text file <ifile>_normalized.dat, with the columns
    1:r 2:rho 3:mass 4:drho/dr."""
    # mass inside each radius: trapezoidal rule for 4*pi*rho*r^2
    y = 4.0 * np.pi * profile * r**2
    mass = np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1])
                                            * np.diff(r))])
    # scales
    R0, M = r[-1], mass[-1]
    rho0 = M / R0**3

    x = r / R0
    rho = profile / rho0
    mass = mass / M
    drhodr = np.gradient(rho, x)  # centered differences; one-sided at ends

    # continue to r = 0 with zero density, with the spacing of the grid
    if x[0] > 0.0:
        n = int(np.ceil(x[0] / (x[1] - x[0])))
        x0 = np.linspace(0.0, x[0], n, endpoint=False)
        zeros = np.zeros(n)
        x = np.concatenate([x0, x])
        rho = np.concatenate([zeros, rho])
        mass = np.concatenate([zeros, mass])
        drhodr = np.concatenate([zeros, drhodr])

    ofile = os.path.splitext(os.path.basename(ifile))[0] + "_normalized.dat"
    header = ("Density profile: from %s\n"
              "scales: R0 = %.12e, M = %.12e, rho0 = M/R0^3 = %.12e\n"
              "1:r[R0]  2:rho[rho0]  3:mass[M]  4:drho/dr[rho0/R0]"
              % (os.path.basename(ifile), R0, M, rho0))
    np.savetxt(ofile, np.column_stack((x, rho, mass, drhodr)), fmt="%.12e",
               header=header)
    print("normalized:  R0 = %.6e, M = %.6e, rho0 = M/R0^3 = %.6e"
          % (R0, M, rho0))
    print("normalized:  written to %s" % ofile)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ifile", help="input HDF5 file")
    parser.add_argument("--r-centering", choices=["nodes", "cells"],
                        default="nodes",
                        help="radial coordinates are grid nodes or cell "
                             "centers [%(default)s]")
    parser.add_argument("--threshold", type=float, default=0.0,
                        help="density above this value counts as non-zero "
                             "for the extents of the matter [%(default)g]")
    parser.add_argument("--plot-profile", action="store_true",
                        help="verify that the density is spherically "
                             "symmetric and plot its radial profile into "
                             "<input>_profile.png")
    parser.add_argument("--output-radial-profile", action="store_true",
                        help="verify that the density is spherically "
                             "symmetric and write its radial profile into "
                             "<input>.dat")
    parser.add_argument("--output-normalized-radial-profile",
                        action="store_true",
                        help="verify that the density is spherically "
                             "symmetric and write its radial profile, "
                             "rescaled to unit radius and unit total mass, "
                             "into <input>_normalized.dat")
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

    if (args.plot_profile or args.output_radial_profile
            or args.output_normalized_radial_profile):
        symmetric = check_symmetry(rho, r, theta, phi)
        # radial profile: average over the solid angle, if not symmetric
        profile = rho[0, 0, :] if symmetric else rho_r / (4.0 * np.pi)
    if args.plot_profile:
        plot_profile(args.ifile, rho, r, symmetric)
    if args.output_radial_profile:
        write_profile(args.ifile, r, profile)
    if args.output_normalized_radial_profile:
        write_normalized_profile(args.ifile, r, profile)


if __name__ == "__main__":
    main()
