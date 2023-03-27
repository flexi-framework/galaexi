## Plane Turbulent Channel Flow
\label{sec:tut_ptcf}

This tutorial describes how to set up and run the Plane-Turbulent-Channel-Flow test case. We will learn how to use the split form DG method to guarantee non-linear stability of the turbulent channel flow. In a second step, we add the sub grid scale model of Smagorinsky combined with Van Driest type damping to run stable wall-bounded turbulent flows with explicit small scale dissipation. The tutorial assumes that you are familiar with the general **FLEXI** and **HOPR** work flow (please finish the previous tutorials first if this sounds strange to you).

### Flow description

The flow is calculated in a plane channel with half-height $\delta=1$, streamwise (x coordinate) length $2\pi$ and span (z coordinate) width $\pi$ with periodic boundaries in the x- and z-directions as well as no-slip walls at the top and the bottom of the domain. As initial conditions an analytical mean turbulent velocity profile is used superimposed with sinus perturbations in the u, v and w velocity components and a constant density of $\rho=1$. The superimposed pertubations lead to rapid production of turbulent flow structures. Since the wall friction would slow down the flow over time, a constant pressure source term imposing a pressure gradient $\frac{dp}{dx}=-1$ is added as a volume source. While the test case is incompressible in principle, we solve it here in a compressible setting. The chosen Mach number with respect to the bulk velocity in the field is $Ma=0.1$ according to the Moser channel test case. In this setting, the wall friction velocity $\tau$ will always be equal to $1$. We can the define a Reynolds number based on the channel half-height and the wall friction velocity as $Re_{\tau}=1/\nu$.

### Compiler options

Make sure that **FLEXI** is compiled with the CMake options listed in the following table.


| Option                          | Value              | Comment      |
| ------------------------------- |:-------------:     | ------------:|
| CMAKE_BUILD_TYPE                | Release            |              |
| FLEXI_EQYNSYSNAME               | navierstokes       |              |
| FLEXI_PARABOLIC                 | ON                 |              |
| LIBS_USE_MPI                    | ON                 |  optional    |
| FLEXI_EDDYVISCOSITY             | ON                 |  optional    |
| FLEXI_NODETYPE                  | GAUSS-LOBATTO      |              |
| FLEXI_SPLIT_DG                  | ON                 |              |
| FLEXI_TESTCASE                  | channel            |              |
| POSTI                           | ON                 |              |
| POSTI_CHANNEL_FFT               | ON                 |              |

Table: CMake options for the plane turbulent channel flow test case simulation. \label{tab:ptcf_cmakeoptions}

For all other CMake options you may keep the default values. Compile the code.

#### Mesh Generation with HOPR

We use a cartesian mesh with 4 cells per direction for the tutorial. The mesh is stretched in the wall-normal direction to accommodate for the straining of the vortexes close to the wall. In case you want to generate other meshes the parameter file for **HOPR** is included in the tutorial directory (*parameter_hopr.ini*),
the default mesh is included. Using 4 cells with a polynomial degree of $N=5$, means we use a large eddy simulation setup of $24$ DOFs per direction.

### Tutorial - Flow at $Re_{\tau}=180$

Copy the ``plane_turbulent_channel_flow`` tutorial folder to your working directory.

        cp -r $FLEXI_TUTORIALS/plane_turbulent_channel_flow .

Step into the folder. In case you do not want to generate the mesh files yourself, a default mesh has already been provided.

#### Preparing the Flow Simulation with FLEXI

The simulation setup is already defined in *parameter_flexi.ini*.

##### Output

In this tutorial we don't look at the flow visualization of the instantaneous state files. Here, we will rather post process consecutive, instantaneous state files with the ``posti_channel_fft`` tool. As an output, we receive mean velocity and Reynolds stress profiles as well as turbulent energy spectra at different locations normal to the channel wall.

##### Interpolation / Discretization parameters

In this tutorial we use the split form DG method to guarantee non-linear stability of the turbulent channel flow simulation. As already specified in the CMake options table \ref{tab:ptcf_cmakeoptions}, the ``FLEXI_SPLIT_DG`` option has to be switched ON in combination with the ``FLEXI_NODETYPE`` ``GAUSS-LOBATTO``. **FLEXI** has several split flux formulations implemented. Therefore, a specific split flux formulation has to be set in the *parameter_flexi.ini* file. In this tutorial the pre-defined split flux formulation by Pirozzoli is used, which results in a kinetic energy preserving DG scheme.

~~~~~~~
! ================================================ !
! SplitDG
! ================================================ !
SplitDG       = PI     ! SplitDG formulation to be used: SD, MO, DU, KG, PI
~~~~~~~

To switch on Smagorinsky's model set the eddyViscType to $1$ in the *paramerter_flexi.ini* file. In addition, the following parameters have to be set. CS is the Smagorinsky constant usually chosen around $0.11$ for wall bounded turbulent flows and the turbulent Prandtl number is commonly set to $0.6$. To ensure the correct behaviour of the eddy viscosity provided by Smagorinsky's model when approaching a wall, Van Driest type damping has to be switched on.

~~~~~~~
! ================================================ !
! LES MODEL
! ================================================ !
eddyViscType = 0       ! Choose LES model, 1:Smagorinsky
VanDriest = T          ! Van Driest damping for LES viscosity (channel flow only)
CS = 0.11              ! Smagorinsky constant
PrSGS = 0.6            ! turbulent Prandtl number
~~~~~~~

#### Running the Simulation and Results

Now run the simulation, either using

~~~~~~~
flexi parameter_flexi.ini
~~~~~~~

or

~~~~~~~
mpirun -np XX flexi parameter_flexi.ini
~~~~~~~

when you want to use more than one processor.
Once the simulation finished state files can be post processed by the ``posti_channel_fft`` tool which was build by the ``POSTI_CHANNEL_FFT`` CMake option. To run the postprocessing, the standard command is

~~~~~~~
posti_channel_fft parameter_channel_fft.ini [State1 State2 ...]
~~~~~~~

where the *parameter_channel_fft.ini* file is given in the tutorial folder and the amount of statefiles is specified by the user. In this tutorial we use all state files with a timestamp between $t=10.0$ and $t=15.0$. As an output you receive three different files. One containing the mean velocity profiles as well as the Reynolds stress profiles and the other two files contain turbulent erergy spectra.
To visualize those files you can run the python script ``plotChannelFFT.py`` in the ``tools/testcases`` folder with the following command in your simulation directory

~~~~~~~
python $FLEXIROOT/tools/testcases/plotChannelFFT.py -p $PROJECTNAME -t $POSTITIME
~~~~~~~

where ``$PROJECTNAME`` specifies the project name specified in the *parameter_flexi.ini* file and ``$POSTITIME`` the timestamp of your output files from the ``posti_channel_fft`` tool.

#### Part I: SplitDG iLES
First, we run **FLEXI** without Smagorinsky's model which we call an implicit LES (iLES), as no explicit sub-grid scale dissipation model is added. The resulting mean velocity and Reynolds stress profiles as well as turbulent energy spectra close to the centre of the channel are given in Figure \ref{fig:Re180_turbulentChannel}.

![Mean velocity and Reynolds stress profiles (left) as well as turbulent energy spectra close to the centre of the channel (right} of an implicit LES at $Re_{\tau}=180$. \label{fig:Re180_turbulentChannel}](tutorials/10_planeTurbulentChannelFlow/Re180_turbulentChannel.png)

#### Part II: SplitDG with explicit LES model
In a second step, we run **FLEXI** with Smagorinsky's model and Van Driest damping which needs to be switched on in the parameter file as described above. The resulting mean velocity and Reynolds stress profiles as well as turbulent energy spectra close to the centre of the channel are given in Figure \ref{fig:Re180_turbulentChannel_Smag}. In comparison to the previous simulation you might recognize the effect of the explicit damping on the Reynolds stress profile $\overline{u'u'}$ close to the maximum, most. To further study the influence of Smagorinsky's model play around with the spatial resolution both in terms of grid resolution as well as the polynomial degree N. You can also increase the Reynolds number to $Re_{\tau}=395$ or $Re_{\tau}=590$ and compare the results to DNS results from Moser et al. [@moser1999direct].

![Mean velocity and Reynolds stress profiles (left) as well as turbulent energy spectra close to the centre of the channel (right} of a LES with Smagorinsky's model and vanDriest damping at $Re_{\tau}=180$. \label{fig:Re180_turbulentChannel_Smag}](tutorials/10_planeTurbulentChannelFlow/Re180_turbulentChannel_Smag.png)

### Performance improvements
\label{sec:tut_ptcf_performance}
FLEXI comes with some advanced optimizations in order to increase its computational efficiency for compute-intesive simulations.
These appear once the flag ``FLEXI_PERFORMANCE=ON`` is set.
The first option is ``FLEXI_PERFORMANCE_OPTLIFT``, which can be activated to optimize the computation of the parabolic terms of the applied equation system.
However, POSTI is not available if this option is enabled.
The option ``FLEXI_PERFORMANCE_PGO=ON`` allows to enable profile-guided optimization (PGO).
For PGO, the executable is first instrumented with some profiling tools by the compiler and then executed on a simple test case.
The generated profiling data can be used by the compiler to identify bottlenecks and hotspots in the code that it cannot spot from the static source code by itself.
The executable is thus compiled a second time using the gathered profiling data to perform these additional optimizations.
In FLEXI, this two-step compilation works as follows.
First, FLEXI is compiled with the following options.

| Option                          | Value              | Comment      |
| ------------------------------- |:-------------:     | ------------:|
| CMAKE_BUILD_TYPE                | Profile/Release    |              |
| FLEXI_EQYNSYSNAME               | navierstokes       |              |
| FLEXI_PARABOLIC                 | ON                 |              |
| LIBS_USE_MPI                    | ON                 |  optional    |
| FLEXI_EDDYVISCOSITY             | ON                 |  optional    |
| FLEXI_NODETYPE                  | GAUSS-LOBATTO      |              |
| FLEXI_SPLIT_DG                  | ON                 |              |
| FLEXI_TESTCASE                  | channel            |              |
| FLEXI_PERFORMANCE               | ON                 |              |
| FLEXI_PERFORMANCE_OPTLIFT       | ON                 |              |
| FLEXI_PERFORMANCE_PGO           | ON                 |              |
| POSTI                           | OFF                |              |

Table: CMake options for the plane turbulent channel flow test case with additional performance flags enabled. \label{tab:ptcf_cmakeoptions_performance}

For the first step, FLEXI havs to be compiled with ``CMAKE_BUILD_TYPE=Profile`` in order to activate the profiling.
Then, FLEXI has to be executed on a simple test case like the freestream tutorial by following the instructions from Section \ref{sec:tut_freestream}.
Finally, FLEXI is compiled a second time, but this time with the build type set to ``CMAKE_BUILD_TYPE=Release`` in order to incorporate the generated profiling data into the compilation process.
Now, FLEXI can be executed as usual and should show a considerable performance improvement in comparison to the previous simulations.
Please be aware, that PGO is currently only supported for the GNU compiler and that this two-step compilation process has to be performed each time the compile options, i.e. the FLEXI executable, are changed.
