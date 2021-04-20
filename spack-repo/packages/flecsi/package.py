# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *
from spack.pkg.lanl_ristra.flecsi import Flecsi

class Flecsi(Flecsi):
    version('flecsph-1', commit='c9718fe208a94adc01771c5df3379ef10bc0f4e2', submodules=True, preferred=False)
    version('flecsph-0', commit='8d4c394a97d0f4eae0a3a772738b454baa892267', submodules=True, preferred=True)

