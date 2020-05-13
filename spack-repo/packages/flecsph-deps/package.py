# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
from spack import *


class FlecsphDeps(Package):

    homepage = "http://flecsi.lanl.com"
    git      = "git@gitlab.lanl.gov:laristra/flecsph.git"

    version('develop', sha256='18d459400558f4ea99527bc9786c033965a3db45bf4c6a32eefdc07aa9e306a6', url="https://www.x.org/archive/individual/util/util-macros-1.19.1.tar.bz2")

    depends_on('cmake@3.12.4:', type='build')
    depends_on('boost@1.70.0: cxxstd=14 +program_options')
    depends_on('mpi')
    depends_on('hdf5@1.8: +mpi')
    depends_on('flecsi backend=mpi')
    depends_on('gsl')
    depends_on('googletest')
    depends_on('ninja')

    def install(self, spec, prefix):
        mod_script = "load_env.sh"

        with open(os.path.join(spec.prefix, mod_script), 'w') as f:
            f.write("# load env")
            f.write("")
            for dep in spec.dependencies(deptype='build'):
                f.write(dep.format())

