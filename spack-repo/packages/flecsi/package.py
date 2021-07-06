from spack import *
# Spack's import hook doesn't support "from spack.pkg.builtin import legion":
from spack.pkg.builtin.flecsi import Flecsi

class Flecsi(Flecsi):
    """
    Additional named versions for flecsi.
    """
    version('2.2', commit='076c39ba276dd162f0a369a8e5d8680e63307a42')
