# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *
from spack.pkg.builtin.flecsi import Flecsi

class Flecsi(Flecsi):

    git = "ssh://git@gitlab.lanl.gov/flecsi/flecsi.git"

    version('2.2.1', commit='fc7640fbc65654c8d65c69aef73a2d5ed7773453', submodules=False, preferred=True)
    version('2.1.0', tag='v2.1.0', submodules=False, preferred=False)


    # All Current FleCSI Releases
    for level in ('low', 'medium', 'high'):
        depends_on('caliper@2.0.1~adiak~libdw', when='@:1.9 caliper_detail=%s' % level)
        depends_on('caliper@2.4.0~libdw', when='@2.0: caliper_detail=%s' % level)

    depends_on('openmpi+legacylaunchers', when='+unit_tests ^openmpi')
    depends_on('legion+shared network=gasnet', when='backend=legion')
    depends_on('legion@cr.9:', when='backend=legion @2.0:')


    # Current FleCSI@:1.9 releases do not support kokkos, omp, or cuda
    conflicts('+kokkos', when='@:1.4')
    conflicts('+openmp', when='@:1.4')
    conflicts('+cuda', when='@:1.4')
    # Unit tests require flog support
    conflicts('+unit_tests', when='~flog')
    # Disallow conduit=none when using legion as a backend
    conflicts('legion conduit=none', when='backend=legion')
    # Due to overhauls of Legion and Gasnet spackages
    #   flecsi@:1.9 can no longer be built with a usable legion
    conflicts('backend=legion', when='@:1.9')

