![logo](doc/flecsph_logo_bg.png)

[![Build Status](https://travis-ci.com/laristra/flecsph.svg?branch=master)](https://travis-ci.com/laristra/flecsph)
[![codecov.io](https://codecov.io/github/laristra/flecsph/coverage.svg?branch=master)](https://codecov.io/github/laristra/flecsph?branch=master)
<!---
[![Quality Gate](https://sonarqube.com/api/badges/gate?key=flecsph%3A%2Fmaster)](https://sonarqube.com/dashboard?id=flecsph%3A%2Fmaster)
--->

# Introduction 

FleCSPH is a multi-physics compact application that exercises FleCSI parallel data structures for tree-based particle methods. In particular, FleCSPH implements a smoothed-particle hydrodynamics (SPH) solver for the solution of Lagrangian problems in astrophysics and cosmology. FleCSPH includes support for gravitational forces using the fast multipole method (FMM).


This project implements smoothed particles hydrodynamics (SPH) method of
simulating fluids and gases using the FleCSI framework.
Currently, particle affinity and gravitation is handled using the parallel
implementation of the octree data structure provided by FleCSI.

We provide several examples of physics problems in 1D, 2D and 3D:

- Sod shock tubes in 1D/2D/3D;
- Noh shock test in 2D/3D;
- Sedov blast waves 2D and 3D;
- airfoil flow in a wind tunnel (2D/3D);
- pressure-induced spherical implosion (2D/3D);
- single and binary stars with Newtonian gravity in 3D.

# Building FleCSPH with Spack

FleCSPH can now be installed as a Spack package. 

In order to install FleCSPH on your machine using spack: 
- Download spack at: github.com/spack/spack 
- Follow installation instructions 
- Use the following command to install core spack utilities:
```{engine=sh}
spack bootstrap
```
- Run:
```{engine=sh}
spack install flecsph 
```
This will build all the dependencies, compile and install FleCSPH. 
In order to use FleCSPH executables simply run: 
```{engine=sh}
spack load flecsph 
```

You will then have access to the generators and the drivers: 
- sodtube\_{1-2-3}d\_generator, sedov\_{1-2-3}d\_generator...
- hydro\_{1-2-3}d, newtonian\_{1-2-3}d...

Sample parameter files and the intial data can be found on the FleCSPH github repository.


## Using Spack in the development workflow (general case)

If you have downloaded FleCSPH from github and working on a development branch, it is very
convenient to use spack to automatically handle the dependencies:

1. Follow the steps above to install FleCSPH with spack. 
This will ensure that all the dependencies are satisfied.
Select your compiler / MPI combination at this step, e.g. use:
```{engine=sh}
spack install flecsph@refactor %gcc@9.1.0 ^openmpi@3.1.4
```
Version `refactor` corresponds to the branch `jloiseau/refactor` on GitLab.
Another version is `master`, it is for the `master` branch on the same repo.

2. To inspect the dependencies:
```{engine=sh}
spack module tcl loads --dependencies flecsph@refactor
```
If this command returns empty, use `spack bootstrap` for tcl.

3. Load the FleCSPH dependencies installed by spack into the ``bash`` environment:
```{engine=sh}
source <(spack module tcl loads --dependencies flecsph@refactor)
```
Unload FleCSPH itself as you will be using your own custom built version:
```{engine=sh}
module unload $(spack module tcl find flecsph@refactor)
```
Inspect your module environment to make sure dependencies have been loaded:
```{engine=sh}
module list
```

4. You can now build your development version with cmake as described below, 
skipping all the dependencies.
cmake should find all the dependencies from what you loaded with spack:
```{engine=sh}
mkdir build; cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=debug \
    -DENABLE_UNIT_TESTS=ON   \
    -DENABLE_DEBUG=OFF       \
    -DLOG_STRIP_LEVEL=1
```

## Precompiled modules on yellow / turquoise clusters

For the new branch (`jloieau/refactor`), there are precompiled dependency modules both
in project directories on turquoise and yellow clusters. You can preload them and skip
compiling the dependencies:

1. Source the file with modules on which FleCSPH depends (compiled with GCC/8.3.0 and OpenMPI/2.1.2):
  
```{engine=sh} 
    # On yellow clusters, snow or grizzly:   
    ssh sn-fey    # or #    ssh gr-fey
    module purge
    source /usr/projects/packages/flecsph/env_gcc-8.3.0_openmpi-2.1.2.sh

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

## On Darwin: using Spack to install FleCSPH

 installation instructions of the jloiseau/refactor branch for Darwin (with GCC/8.2.0 and MPICH/3.2.1).

1. Clone spack and run bootstrap:
 
```{engine=sh}                
    cd ~/src
    git clone --recursive git@github.com:spack/spack
    source $HOME/src/spack/share/spack/setup-env.sh # add this to ~/.bashrc
    spack bootstrap
```

2. Add the custom spack-repo for FleCSPH:

```{engine=sh}                
    git clone git@gitlab.lanl.gov:laristra/flecsph
    git checkout jloiseau/refactor
    spack repo add ~/num/FleCSPH/flecsph/spack-repo
```

3. Load compiler and cmake modules:

```{engine=sh}                
    module load cmake/3.12.4 gcc/8.2.0
```

4. Create file `~/.spack/linux/packages.yaml` with the following content:

```{engine=sh}                
-- >>> -----------------------------------------------
packages:
    cmake:
        modules:
            cmake@3.12.4: cmake/3.12.4
    mpich:
        modules:
            mpich@3.2.1-gcc_8.2.0: mpich-slurm/3.2.1-gcc_8.2.0
    all:
        compiler: [gcc@8.2.0]
        providers:
            cmake: [cmake@3.12.4]
            mpi: [mpich@3.2.1-gcc_8.2.0]
-- <<< -----------------------------------------------
```

5. Install FleCSPH:

```{engine=sh}                
    spack install flecsph@refactor %gcc@8.2.0 ^mpich@3.2.1
    # repeat if fails
```

6. Setup the environment:

```{engine=sh}                
    module purge
    module load cmake/3.12.4 gcc/8.2.0
    source <(spack module tcl loads --dependencies flecsph@refactor)
    module unload $(spack module tcl find flecsph)
```

7. Compile and test:

```{engine=sh}                
    cmake .. \
         -DCMAKE_BUILD_TYPE=debug \
         -DENABLE_UNIT_TESTS=ON   \   
         -DENABLE_DEBUG=OFF       \
         -DLOG_STRIP_LEVEL=1
    make -j
    make test
```



# Building FleCSPH manually

FleCSPH can be installed anywhere in your system; to be particular, below we
assume that all repositories are downloaded in FLECSPH root directory `${HOME}/FLECSPH`.

## Suggested directory structure

We recommend to use an isolated installation of FleCSPH and FleCSI, such that the software and all their
dependencies in a separate directory, with the following directory structure:

```{engine=sh}
  ${HOME}/FLECSPH
  ├── flecsi
  │   └── build
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
- GSL library 

## FleCSI

Clone FleCSI repo at the `master` branch (default).
Checkout submodules recursively, then configure as below:

```{engine=sh}    
   export CMAKE_PREFIX_PATH=${HOME}/FLECSPH/local
   cd $HOME/FLECSPH
   git clone --recursive git@github.com:laristra/flecsi.git
   cd flecsi
   git submodule update --recursive
   mkdir build ; cd build
   cmake .. \
       -DCMAKE_INSTALL_PREFIX=$CMAKE_PREFIX_PATH  \
       -DENABLE_MPI=ON                            \
       -DENABLE_MPI_CXX_BINDINGS=ON               \
       -DENABLE_OPENMP=ON                         \
       -DCXX_CONFORMANCE_STANDARD=c++17           \
       -DENABLE_LOG=ON                           \
       -DFLECSI_RUNTIME_MODEL=mpi                 \
       -DENABLE_FLECSIT=OFF                       \
       -DENABLE_FLECSI_TUTORIAL=OFF               
```    

In this configuration, MPI is used as FleCSI backend.
If you want to use other FleCSI backends (Legion, HPX), you will need to install them separately: see https://github.com/laristra/flecsi-third-party for further info.

In a final step, build and install:

    % make -j install

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

A new evolution driver must have a `main.cc` and `main_driver.cc` files. Do not edit
`main.cc`, because FleCSI expects certain format of this file. It is easier to start
by copying existing files to your folder under `app/drivers`. Include cmake
targets with different dimensions using examples in `app/drivers/CMakeLists.txt`.

Make sure to document your subproject in a corresponding `README.md` file
that describes the problem you want to run. In order to get all files easily and
correctly, you can copy them from other subprojects such as `sodtube` or `hydro`.

# For developers
Please refer to the following page:
[Development Guidelines](https://github.com/laristra/flecsph/blob/doc/hlim/doc/development.md)

## Style guide

FleCSPH follows the FleCSI coding style, which in turn follows (in general) the Google coding conventions.
FleCSI coding style is documented here:
https://github.com/laristra/flecsi/blob/master/flecsi/style.md

# Logs

Cinch Log is the logging tool for this project.
In order to display log set the environment variable as:
```bash
export LOG_ENABLE_STDLOG=1
```

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
