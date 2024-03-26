[![logo](https://numericsresearchgroup.org/images/icons/galexi.svg "Galexi")][flexi]


[![license](https://img.shields.io/github/license/flexi-framework/flexi.svg?maxAge=2592000 "GPL-3.0 License")](LICENSE.md)
[![doi](https://img.shields.io/badge/DOI-10.1016/j.camwa.2020.05.004-blue "DOI")](https://doi.org/10.1016/j.camwa.2020.05.004)
[![youtube](https://img.shields.io/badge/YouTube-red?logo=youtube "YouTube")](https://www.youtube.com/@nrgiag8633)
[![userguide](https://img.shields.io/badge/Userguide-silver "Userguide")][userguide]
[![gallery](https://img.shields.io/badge/Gallery-teal "Gallery")][gallery]

# About

[GALÆXI][flexi] is a GPGPU-enabled extension of FLEXI, a high-order numerical framework for solving PDEs, with a special focus on Computational Fluid Dynamics.
[GALÆXI][flexi] is based on the Discontinuous Galerkin Spectral Element Method (DGSEM), which allows for high-order of accuracy 
and fully unstructured hexahedral meshes. The solver is parallelized very efficiently for large-scale applications using MPI-aware CUDA Fortran and has demonstrated excellent scaling on over 1000 NVIDIA GPUs. Moreover, [GALÆXI][flexi] comes with a capable pre- and postprocessing suite that enables complex
simulation setups up to the finished visualization.

[GALÆXI][flexi] has been developed by the [Numerics Research Group (NRG)][nrg] founded by Prof. Claus-Dieter Munz and currently
lead by Prof. Andrea Beck at the Institute of Aerodynamics and Gasdynamics at the University of Stuttgart, Germany.

You can find detailed installation instructions, extensive documentation and
several tutorial cases for GALÆXI [here][flexi].

GALÆXI is Copyright (C) 2016, Prof. Claus-Dieter Munz and is released under the **GNU General Public License v3.0**.
For the full license terms see the included [license file](LICENSE.md).

Numerous people have worked on and with GALÆXI over the last years.
We would like to thank all these [contributors](CONTRIBUTORS.md) for their efforts they spent on building GALÆXI.
 
In case you have questions regarding GALÆXI or want to contribute yourself
by either reporting bugs, requesting features or adding somthing
different to the project, feel free to open an issue or pull request.

# Cite
GALÆXI is a scientific project. If you use GALÆXI for publications or
presentations in science, please support the project by citing it.
As general reference, please cite
```
Krais, N., Beck, A., Bolemann, T., Frank, H., Flad, D., Gassner, G., Hindenlang, F., Hoffmann, M., Kuhn, T., Sonntag, M., Munz, C.-D.
FLEXI: A high order discontinuous Galerkin framework for hyperbolic–parabolic conservation laws,
Computers & Mathematics with Applications, 81, (2021) 186-219
```
or use the following Bibtex files

    @article{flexi,
      title = {{FLEXI}: {A} high order discontinuous {G}alerkin framework for hyperbolic-parabolic conservation laws},
      journal = {Computers \& Mathematics with Applications},
      volume = {81},
      pages = {186-219},
      year = {2021},
      doi = {https://doi.org/10.1016/j.camwa.2020.05.004},
      author = {Nico Krais and Andrea Beck and Thomas Bolemann and Hannes Frank and David Flad and Gregor Gassner and Florian Hindenlang and Malte Hoffmann and Thomas Kuhn and Matthias Sonntag and Claus-Dieter Munz},
    }

To refer to specific applications and features, you can also cite the appropriate paper from this [list][publications].

# Quick Start Guide
For a more detailed installation instructions, please see the documention [here][userguide].

GALÆXI is tested for various Linux distributions including Ubuntu, OpenSUSE, CentOS and Arch. It currenly only supports offloading to NVIDIA GPUs. Also, compilation is currently only possible with the NVIDIA HPC SDK compilers (nvfortran).
For installation you require the following dependencies:

| Package          | Required | Installed by GALÆXI|
|:-----------------|:--------:|:------------------:|
| Git              |      x   |                    |
| CMake            |      x   |                    |
| C/C++ Compiler   |      x   |                    |
| Fortran Compiler |      x   |                    |
| LAPACK           |      x   |      x             |
| HDF5             |      x   |      x             |
| MPI              |     (x)  |                    |
| CUDA             |      x   |                    |

The MPI library is only required for running parallel simulations on multiple ranks and the HDF5 and LAPACK libraries
can be installed automatically during the GALÆXI build process.
The names of the packages and the package manager might differ depending on the specific distribution used.

### Getting the code
Open a terminal, download GALÆXI via git and optionally export the GALÆXI directory:

    git clone https://github.com/flexi-framework/galaexi.git
    export GALAEXI_DIR="$(pwd)/galaexi"

### Compiling the code
Enter the GALÆXI directory, create a build directory and use CMake to configure and compile the code

    cd $GALAEXI_DIR
    mkdir build; cd build
    cmake ../
    make

The executable `galaexi` is now contained in the GALÆXI directory in `build/bin/`.
Custom configurations of the compiler options, dependencies and code features can be set using

    ccmake ../

### Running the code
Navigate to the directory of the tutorial **naca0012** and run GALÆXI

    cd $GALAEXI_DIR/tutorials/naca0012
    $GALAEXI_DIR/build/bin/galaexi parameter_flexi_navierstokes.ini

# Used libraries
GALÆXI uses several external libraries as well as auxiliary functions from open source projects, including:
* [HDF5](https://www.hdfgroup.org/)
* [MPI](https://www.mcs.anl.gov/research/projects/mpi/)
* [LAPACK](https://www.netlib.org/lapack/)
* [OpenMP](https://www.openmp.org/)
* [FFTW](https://www.fftw.org/)
* [CMake](https://cmake.org/)
* [Reggie2.0](https://github.com/reggie-framework/reggie2.0/)
* [PAPI](https://icl.cs.utk.edu/papi/)
* [CUDA](https://developer.nvidia.com/cuda-toolkit)

[nrg]:           https://numericsresearchgroup.org/index.html
[flexi]:         https://numericsresearchgroup.org/flexi_index.html
[publications]:  https://numericsresearchgroup.org/publications.html#services
[userguide]:     https://numericsresearchgroup.org/userguide/userguide.pdf
[gallery]:       https://numericsresearchgroup.org/gallery.html#portfolio
[youtube]:       https://www.youtube.com/@nrgiag8633 
