![logo](doc/flecsph_logo_bg.png)

[![Build Status](https://travis-ci.com/laristra/flecsph.svg?branch=master)](https://travis-ci.com/laristra/flecsph)
[![codecov.io](https://codecov.io/github/laristra/flecsph/coverage.svg?branch=master)](https://codecov.io/github/laristra/flecsph?branch=master)
<!---
[![Quality Gate](https://sonarqube.com/api/badges/gate?key=flecsph%3A%2Fmaster)](https://sonarqube.com/dashboard?id=flecsph%3A%2Fmaster)
--->

# Introduction 

FleCSPH is a multi-physics compact application for tree-based particle methods. In particular, FleCSPH implements a smoothed-particle hydrodynamics (SPH) solver for the solution of Lagrangian problems in astrophysics and cosmology. FleCSPH includes support for gravitational forces using the fast multipole method (FMM). It uses MPI for distributed-memory parallelism and an internal octree data structure for neighbor finding and gravity.


This project implements smoothed particles hydrodynamics (SPH) method of
simulating fluids and gases.
Currently, particle affinity and gravitation is handled using the parallel
implementation of the octree data structure.

We provide several examples of physics problems in 1D, 2D and 3D:

- Sod shock tubes;
- Noh implosion test;
- Sedov blast wave;
- single and binary stars with Newtonian gravity in 3D;
- Kelvin-Helmholtz instability setup;
- Rayleigh-Taylor instability setup;
- airfoil flow in a wind tunnel.

# Building FleCSPH with Spack

FleCSPH can now be installed as a Spack package. 

In order to install FleCSPH on your machine using spack, follow these steps.
These instructions work for Darwin -- please modify them as necessary for your machine.

1. Clone the spack git repository and load spack environment:
```sh
  git clone https://github.com/spack/spack.git
  source $HOME/src/spack/share/spack/setup-env.sh 
```

2. Use the following command to install the Linux `module` utility:
```sh
  spack install environment-modules
```

3. Clone FleCSPH; purge modules and load `cmake` and GCC compilers (>8.2 version is required).
```sh
  git clone git@gitlab.lanl.gov:laristra/flecsph
  cd flecsph
  module purge
  module load cmake/3.17.0 gcc/9.3.0
```

4. Modify `spack` configuration files, `~/.spack/linux/compilers.yaml` and 
`~/.spack/linux/packages.yaml` to point to the correct version of `cmake`, 
`gcc` and MPI (`openmpi` on Darwin is preferred).

File `~/.spack/linux/compilers.yaml`:
```yaml
compilers:
- compiler:
    spec: gcc@9.3.0
    paths:
      cc: /projects/opt/x86_64/gcc/9.3.0/bin/gcc
      cxx: /projects/opt/x86_64/gcc/9.3.0/bin/g++
      f77: /projects/opt/x86_64/gcc/9.3.0/bin/gfortran
      fc: /projects/opt/x86_64/gcc/9.3.0/bin/gfortran
    flags: {}
    operating_system: centos7
    target: x86_64
    modules: []
    environment: {}
    extra_rpaths: []
```

File `~/.spack/linux/packages.yaml`:
```yaml
packages:
    cmake:
        modules:
            cmake@3.17.0: cmake/3.17.0
    mpich:
        modules:
            mpich@3.2.1-gcc_8.2.0: mpich-slurm/3.2.1-gcc_8.2.0
    openmpi:
        modules:
            openmpi@4.0.3: openmpi/4.0.3-gcc_9.3.0
    all:
        compiler: [gcc@9.3.0]
        providers:
            cmake: [cmake@3.17.0]
            mpi: [openmpi@4.0.3]
```

5. Finally, run FleCSPH installation:
```sh
spack install flecsph %gcc@9.3.0 ^openmpi@4.0.3
```
This will take advantage of the existing modules, build all the missing dependencies (such as `google-test`), compile and install FleCSPH. 

6. In order to use FleCSPH executables, simply run: 
```{engine=sh}
spack load flecsph 
```
You will then have access to the generators and the drivers: 
- sodtube\_[123]d\_generator, sedov\_[123]d\_generator...
- hydro\_[123]d, newtonian\_3d, wvt\_[123]d...

Sample parameter files and the intial data can be found in the `data` subdirectory.


## Installation with Spack in the development workflow

If you have downloaded FleCSPH from github and working on a development branch, it is
convenient to use spack to automatically handle the dependencies:

1. -5. Follow the steps 1-5 above to install FleCSPH with spack. 
This will ensure that all the dependencies are satisfied.

6. To inspect the dependencies:
```sh
spack module tcl loads --dependencies flecsph
```

7. Load the FleCSPH dependencies installed by spack into the ``bash`` environment:
```{engine=sh}
source <(spack module tcl loads --dependencies flecsph)
```
Unload FleCSPH itself as you will be using your own custom built version:
```{engine=sh}
module unload $(spack module tcl find flecsph)
```
Inspect your module environment to make sure dependencies have been loaded:
```{engine=sh}
module list
```

8. You can now build your development version with cmake as described below, 
skipping all the dependencies.
cmake should find all the dependencies from what you loaded with spack:
```{engine=sh}
mkdir build; cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=debug -DENABLE_UNIT_TESTS=ON   \
    -DENABLE_DEBUG=OFF       -DLOG_STRIP_LEVEL=1
make -j
make test
```

9. It is convenient to create a shell script which loads all precompiled 
modules before working on FleCSPH development.
```sh
echo "# Loads FleCSPH environment" > load_flecsph_env.sh
echo "module purge"                           >> load_flecsph_env.sh
echo "module load cmake/3.17.0 gcc/9.3.0"     >> load_flecsph_env.sh
spack module tcl loads --dependencies flecsph >> load_flecsph_env.sh
```
Then, you can run `source load_flecsph_env.sh` to load the dependencies.

## Precompiled modules on yellow / turquoise clusters

There are precompiled dependency modules both in project directories on 
turquoise and yellow clusters. You can preload them and skip compiling 
the dependencies:

1. Source the file with modules on which FleCSPH depends (compiled with GCC/8.3.0 and OpenMPI/2.1.2):
  
```{engine=sh} 
    # On yellow clusters, snow or grizzly:   
    ssh sn-fey    # or #    ssh gr-fey
    module purge
    source /usr/projects/packages/flecsph/snow_env_gcc-9.2.0_openmpi-3.1.5

    # On turquoise clusters badger or grizzly:
    ssh ba-fe     # or #    ssh gr-fe
    module purge
    source /usr/projects/packages/flecsph/env_gcc-8.3.0_openmpi-2.1.2.sh
``` 
         
2. Change to your FleCSPH development directory and start configuring:
```{engine=sh}          
    cd flecsph/build
    ccmake .. -DCMAKE_BUILD_TYPE=debug -DENABLE_UNIT_TESTS=ON -DENABLE_DEBUG=OFF -DLOG_STRIP_LEVEL=1
```
               
3. In ccmake interface, press 't' and correct the flags `MPI_CXX_LINK_FLAGS` and `MPI_C_LINK_FLAGS`: 
   replace
```{engine=sh}                
   -Wl,-rpath -Wl,/usr/lib64
```                      
   with
```{engine=sh}                
   -Wl,-rpath -Wl,/usr/projects/hpcsoft/toss3/common/x86_64/gcc/8.3.0/lib64 
```
   Save and exit.

4. Compile and test your build:
```{engine=sh}
   make -j
   make test
```

## Spack mirrors in the project space on yellow and turquoise

If you want to compile your own Spack modules on Yellow or Turquoise, where 
access to some Internet repositories is restricted, you can use these shared 
mirrors:
```{engine=sh}
 - on yellow: /usr/projects/packages/flecsph/spack_mirror
 - on turquoise: /turquoise/usr/projects/nsmergers/spack/mirror
```

To use mirrors, login to turquoise or yellow and add the mirrors:
 - on yellow: `spack mirror add yellow_mirror file:///usr/projects/packages/flecsph/spack_mirror`
 - on turquoise: `spack mirror add turq_mirror file:///turquoise/usr/projects/nsmergers/spack/mirror`

Mirror directories contain tarballs for various packages. With mirrors added, Spack will resort 
to those tarballs if it cannot reach their standard location on the Internet.
Mirrors are open for writing within the group 'nsmergers'. If some packages are missing, 
you can copy them to the mirrors as described 
[here](https://spack.readthedocs.io/en/latest/mirrors.html).



# Building FleCSPH manually

FleCSPH can be installed anywhere in your system; to be particular, below we
assume that all repositories are downloaded in FLECSPH root directory `${HOME}/FLECSPH`.

## Suggested directory structure

We recommend to use an isolated installation of FleCSPH with the software and all its
dependencies in a separate directory, with the following directory structure:

```{engine=sh}
  ${HOME}/FLECSPH
  ├── flecsph
  │   ├── build
  │   └── third-party-libraries
  └── local
      ├── bin
      ├── include
      ├── lib
      ├── lib64
      └── share
```

Below we use `${HOME}/FLECSPH/local` for an installation directory.
Make sure to set your CMAKE prefix to this location:

    % export CMAKE_PREFIX_PATH=${HOME}/FLECSPH/local

## Prerequisites

You will need the following tools:

- C++17 - capable compiler, such as gcc version >= 7;
- git version > 2.14;
- MPI libraries, compiled with the gcc compiler above and multithread support
  (`--enable-mpi-thread-multiple` for OpenMPI and
   `--enable-threads=multiple` for MPICH);
- cmake version > 3.7;
- boost library version > 1.59;
- Python version > 2.7.
- HDF5 compiled with parallel flag version > 1.8

## FleCSPH

Clone the master branch from the FleCSPH git repo:
```{engine=sh}
   cd $HOME/FLECSPH
   git clone --recursive git@github.com:laristra/flecsph.git
```    

### Building FleCSPH

Configure command:

```{engine=sh}
   # in ${HOME}/FLECSPH/build:
   export CMAKE_PREFIX_PATH=${HOME}/FLECSPH/local
   cmake .. \
       -DCMAKE_BUILD_TYPE=debug \
       -DENABLE_UNIT_TESTS=ON   \
       -DENABLE_DEBUG=OFF       \
       -DLOG_STRIP_LEVEL=1
```

Build and install:

    % make -j install


### Building FleCSPH on various architectures
Architecture-/machine-specific notes for building FleCSPH are collected in
[doc/machines](https://github.com/laristra/flecsph/tree/master/doc/machines).
If you succeeded in compiling and running FleCSPH on new architectures,
please do not hesitate to share your recipe.
We appreciate user contributions.

# Running FleCSPH applications
Current FleCSPH contains several initial data generators and two evolution
drivers: `hydro` and `newtonian`. Initial data generators are located in
`app/id_generators/`:

- `sodtube`: 1D/2D/3D sodtube shock test;
- `sedov`: 2D and 3D Sedov blast wave;
- `noh`: 2D and 3D Noh implosion test.
- etc.

Evolution drivers are located in `app/drivers`:

- `hydro`: 1D/2D/3D hydro evolution without gravity;
- `newtonian`: 3D hydro evolution with self-gravity.

To run a test, you also need an input parameter file, specifying parameters of the
problem. Parameter files are located in `data/` subdirectory. Running an
application consists of two steps:

- generating initial data;
- running evolution code.

For instance, to run a `sodtube` 1D shock test, do the following (assuming
you are in your build directory after having successfully built FleCSPH):
```{engine=sh}
  cp ../data/sodtube_t1_n1000.par sodtube.par
  # edit the file sodtube.par to adjust the number of particles etc.
  app/id_generators/sodtube_1d_generator sodtube.par
  app/driver/hydro_1d sodtube.par
```

# Creating your own initial data or drivers
You can add your own initial data generator or a new evolution module under
`app/id_generators` or `app/drivers` directories. Create a directory with
unique name for your project and modify `CMakeLists.txt` to inform the cmake
system that your project needs to be built.

A new initial data generator usually has a single `main.cc` file and an optional
include file. You can use existing interfaces for lattice generators or equations
of state in the `include/` directory.
The file `app/drivers/include/user.h` defines the dimensions of your problem, both
for initial data generators and for the evolution drivers.
This is done via a compile-time macro `EXT_GDIMENSION`, which allows users to have
the same source code for different problem dimensions. Actual dimension is set at
compile time via the `target_compile_definitions` directive of cmake, e.g.:
```
   target_compile_definitions(sodtube_1d_generator PUBLIC -DEXT_GDIMENSION=1)
   target_compile_definitions(sodtube_2d_generator PUBLIC -DEXT_GDIMENSION=2)
```

A new evolution driver must have a `main_driver.cc` file. It is easier to start
by copying existing files to your folder under `app/drivers`. Include cmake
targets with different dimensions using examples in `app/drivers/CMakeLists.txt`.

Make sure to document your subproject in a corresponding `README.md` file
that describes the problem you want to run. In order to get all files easily and
correctly, you can copy them from other subprojects such as `sodtube` or `hydro`.

# For developers
Please refer to the following page:
[Development Guidelines](https://github.com/laristra/flecsph/blob/master/doc/development.md)

# Logs

In the code, you can set the level of output from trace(0) and info(1) to warn(2), error(3) and fatal(4).
You can then control the level of output at compile time by setting the flag `LOG_STRIP_LEVEL`:
by default it is set to 0 (trace), but for simulations it is perhaps preferrable to set it to 1 (info).
```cpp
log_one(trace) << "This is verbose output  (level 0)" << std::endl;
log_one(info) << "This is essential output (level 1)" << std::endl;
log_one(warn) << "This is a warning output (level 2)" << std::endl;
log_one(fatal) << "Farewell!" << std::endl;
```

For further details, refer to the documentation at:
https://github.com/laristra/cinch/blob/master/logging/README.md

# Contacts

If you have any questions or concerns regarding FleCSPH, please contact Julien Loiseau (jloiseau@lanl.gov), 
Oleg Korobkin (korobkin@lanl.gov) and/or Hyun Lim (hyunlim@lanl.gov)
