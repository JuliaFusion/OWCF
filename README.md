![The logo of the orbit weight computational framework](misc/logo.gif)

## 0. This is the Orbit Weight Computational Framework (OWCF)

It enables effective and efficient computation of weight functions using the Julia programming language (coupled to other programming languages, e.g. Python). As of the current version (v1.7.0), the OWCF can compute the following types of weight functions:
- Velocity-space weight functions (2D). Computed by feeding (E,p) points into a synthetic diagnostic code of your choice, while keeping (R,z) constant. It can output and save the 2D weight functions in both (E,p) and (vpara,vperp)
- Orbit-space weight functions (3D). Computed using the guiding-centre computational framework (GCCF)(written in the Julia language) developed by Luke Stagner with additions by Henrik Järleblad. The GCCF is comprised of the Julia packages EFIT.jl, Equilibrum.jl, GuidingCenterOrbits.jl and OrbitTomography.jl. The OWCF computes fast-ion orbits in the (E,pm,Rm) orbit space, and feeds them as weighted (E,p,R,Z) data into the forward model of your choice
- (E,p,R,z) weight functions (4D). Computed by feeding (E,p,R,z) points into a synthetic diagnostic code of your choice. It can output and save the 4D weight functions in both (E,p,R,z) and (vpara,vperp,R,z) coordinates
 
The OWCF also has the ability to compute expected diagnostic spectra, given a fusion reaction, pertaining fast and bulk (thermal) reactant distribution and a toroidally symmetric magnetic equilibrium. This can be done via the script calcSpec.jl. In fast-ion tomography, such a spectrum of (synthetic) measurements is identified as S. In fast-ion tomography, the main equation is

> `S = WF`

where S is a vector with the measurements (synthetic or experimental), W is the weight function that captures all of the necessary physics in a simple matrix, and F is the fast (ion) reactant distribution reshaped as a vector. To clarify, the OWCF thus has the ability to compute both S and W (both not F, i.e. the OWCF has no Fokker-Planck equation solver). F needs to be provided as e.g. a TRANSP distribution. You can also create your own fast-ion distribution in 2D, 3D or 4D coordinates. This is done with the OWCF/extra/createCustomFIDistr.jl tool. As long as you save the fast-ion distribution data in a file with the data keys explained in e.g. ps2os.jl, it will be usable as F input.

The OWCF_manual.pdf contains more detailed (but possibly outdated) information than this simple README. Please consult the OWCF_manual.pdf document for more information on the OWCF. More information about the OWCF can also be obtained by reading the following publication:

> H. Järleblad *et al*, *A framework for synthetic diagnostics using energetic-particle orbits in tokamaks*, Comp. Phys. Comm. **294**, 108930 (2024)  

## 1. To install
### 1.1 Julia
To install Julia, we recommend that you install the JuliaUp version manager, to be able to easily switch between, and manage, different Julia versions when needed. To download and install JuliaUp on a Unix-type system (Linux or Mac), copy and run the following command in your terminal:

> `curl -fsSL https://install.julialang.org | sh`

To download and install Julia on a Windows-type system, copy and run the following command in your command prompt (e.g. Powershell):

> `winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore`

After downloading and installing JuliaUp, to ensure that the OWCF will be working correctly **you should install Julia version 1.11.7 or 1.11.6**. To do this, run the following commands in your command-line terminal:

> `juliaup add 1.11.7`  
> `juliaup default 1.11.7`

Afterwards, to start Julia, simply type `julia` in your command-line terminal.

### 1.2 Python
Since the OWCF can utilize the DRESS code [J. Eriksson *et al* 2016 Comp. Phys. Comm. **199** 40-46], written in Python, when computing expected neutron and gamma-ray product energies for fusion reactions, a Python executable needs to be provided to connect Julia with Python. In addition, to ensure correct functionality, **we recommend Python version 3.12.10**.

### 1.3 Python packages
The **numpy**, **scipy** and **pandas** Python packages need to be installed for your Python executable. Furthermore, if you want to use TRANSP [Breslau J., Gorelenkova M., Poli F., Sachdev J., Pankin A. and Perumpilly G. 2018 TRANSP (https://doi.org/10.11578/dc.20180627.4)] output files in .cdf file format as input files to load thermal plasma data from when computing e.g. neutron/gamma-ray energies, the **netCDF4** package also needs to be installed for your Python executable. To install these Python packages, run the following commands in your command-line terminal (replace 'python' and 'pip' with your specific Python executable and Python package manager, if you have any):

> `python -m pip install numpy`  
> `python -m pip install scipy`  
> `python -m pip install pandas`  
> `python -m pip install netCDF4`  

### 1.4 OWCF
To install the OWCF, open a command-line terminal, navigate to the OWCF folder with the `cd` command and type 

> `julia install_OWCF.jl`

and follow the instructions. **PLEASE NOTE!** During the OWCF installation, you need to specify the same Python executable as you used in section 1.3 (see above).

### 1.5 Jupyter notebook
If you would like to use the interactive OWCF apps to analyze and visualize data, you have to install Jupyter notebook. This can be done by running the following command in your command-line terminal (replace 'python' and 'pip' with your specific Python executable and Python package manager, if you have any):

> `python -m pip install notebook`

The Python command ('python') should match the Python command that you specified in section 1.3 and 1.4.

> #### PLEASE NOTE!
> For some scripts, the OWCF currently relies on the [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) package. This package will **not** function correctly if the Python environment that you want to connect is a virtual environment created by `conda`. Therefore, as of the current version, if you are using the OWCF together with a Python environment created by `conda`, correct functionality **cannot** be guaranteed.

## 2. Prior to first use
Before using the OWCF for the first time after downloading and installing, make sure to run the `OWCF/tests/run_tests.jl` script. This script will run tests to check that the OWCF was installed correctly and can be expected to function correctly. Prior to running, please edit the `run_tests.jl` script and modify the input arguments at the top as you see fit. Then do the following:

> *navigate to the OWCF folder using a command-line terminal*  
> `cd tests`  
> `julia run_tests.jl`  

If the `run_tests.jl` script prints `0` errors upon completion, you are good to go!

## 3. How to use
To use the OWCF, find out which script you need. Copy the template start file from the templates/ folder into the main OWCF/ folder. Specify the required inputs. Execute. Get your output results file. Visualize the results using any of the OWCF apps in apps/ or via a plot function in extra/gui.jl Or manually, by writing the necessary code in e.g. scratch.jl.

Some scripts are executed without a template file. The input variables are then specified by modifying the script directly.

**If executing on a SLURM computational cluster, copy the .sh submit template you need from templates/ into the main OWCF/ folder. Modify it as needed. Submit it to the workload manager.**

## 4. Currently supported forward models
A 'forward model' is technical jargon and refers to the process of obtaining synthetic measurements given input distributions in a phase-space of your choice.
For fusion reactions, the distributions are usually a fast (ion) reactant distribution and a bulk (thermal) reactant distribution.

- The DRESS code [J. Eriksson *et al* 2016 Comp. Phys. Comm. **199** 40-46]  
The DRESS code is a code for computing the energy and momentum of fusion reaction products, given the collision and fusing of two reacting nuclei. It also has the possibility to load, handle and use TRANSP data. Given a fusion reaction, it computes the resulting energy spectrum of the emitted particle (e.g. neutrons, gammas). It incorporates diagnostic lines-of-sight via the input of data from the LINE21 code. The DRESS code is modular in the sense that any two-body fusion reaction can be added to the code, simply via defining the corresponding cross-sectional data. As of this version of the OWCF and DRESS code, the DRESS code can compute one-step and two-step fusion reactions.
- Projected velocities.
The projected velocities of the (fast) ion can be used to compute weight functions that quantify the sensitivity of (fast-ion) diagnostics in phase space. The 
velocity vectors of the (fast) ion are projected onto the direction running parallel to the diagnostic LOS. The projected velocity is the underlying quantity 
governing the measurements in several (fast-ion) diagnostics such as collective Thomson scattering, neutron emission spectroscopy, gamma-ray spectroscopy etc.
- Analytic beam-target equations [A. Valentini *et al* 2025 Nucl. Fusion **65** 026001, A. Valentini *et al* 2014 Rev. Sci. Instrum. **95** 083551, A. Valentini *et al* 2025 Nucl. Fusion **65** 046031]
If the temperature of one of the reactants in a fusion reaction can be approximated as zero, analytic equations may be used to compute synthetic diagnostic spectra (and by extension weight functions). The analytic equations are based on the conservation of energy and momentum. Relativistic speeds are supported. Since no Monte Carlo (MC) computations are needed, this forward model is (in theory) orders of magnitude faster than corresponding MC models.

## 5. Currently supported magnetic equilibrium files

- .eqdsk/.geqdsk files. These kind of files are used as input files to TRANSP. They contain information about the magnetic flux function and the tokamak geometry (first wall, divertor etc). For a full breakdown of what an .eqdsk/.geqdsk file usually contains, and its structure, please examine the OWCF/misc/eqdsk_file_breakdown.pdf file.  

- .jld2 file. If the magnetic equilibrium is specified with a .jld2 file, as of the current version of the OWCF (v1.6.0), it has to be an output file of the OWCF/extra/createCustomMagneticEquilibrium.jl script. Or have an equivalent structure and data content. The createCustomMagneticEquilibrium.jl script provides the possibility to create customized Solov'ev magnetic equilibria.

## 6. TRANSP .cdf/.CDF files

- TRANSP shot files. They are usually saved on the form XXXXXAXX.cdf where 'X' can be any (different) numbers and 'A' any letter. E.g. 94701V01.cdf, or 99965K71.cdf. These files is essentially the output file of the TRANSP run with the same ID (e.g. 94701V01). They contain a multitude of data (hence, they are usually quite large) and the OWCF only needs them for their thermal plasma density and temperature profiles.  
- TRANSP fast-ion shot files. They are usually computed using the NUBEAEM module of TRANSP. These files contain the expected fast-ion distributions for the time points of interest. They contain fast-ion data but also time window information data that the OWCF uses to select the corresponding thermal plasma density and temperature profiles. These kind of files can be provided as input to the OWCF/extra/getEpRzFIdistrFromTRANSP.jl script, which will convert it to a .jld2 file with the fast-ion specified on a regular (E,p,R,z) grid.

## 7. Currently supported options for modelling diagnostic sightlines
A sightline is technical jargon for the voxels observed by a diagnostic detector. Synonyms include viewing cone and line-of-sight (LOS).
The OWCF currently supported the following options for LOS:

- LINE21 code (might not be actual name of code) output files. These have diagnostic sightline data 
organized into columns. The columns are, in order: x, y, z, c (angular weight), volume, ux (x component of vector
pointing towards the detector at the given x,y,z position), uy, uz, R, phi, solid angle. Talk to Massimo Nocente (massimo.nocente@mib.infn.it)
or Jacob Eriksson (jacob.eriksson@physics.uu.se) to acquire LINE21 files for diagnostic viewing cones.  
- OWCF custom LOS. The user can use the OWCF/extra/createCustomLOS.jl script to create a custom LOS of their own. This can be done in several different ways and the reader is referred to the description at the start of the createCustomLOS.jl script for a detailed explanation. The output file will be a file identical in structure to a LINE21 output file. The file extension will be .vc but the file can be read like a .txt file.
- No viewing cone. If no viewing cone file is specified, sperical emission is assumed (4*pi s.r.)

## 8. Overview of OWCF main scripts
### 8.1 The main tools (scripts) of the OWCF include:

- calc2DWeights.jl. Compute two-dimensional velocity-space weight functions for an (E,p) grid, at a given (R,z) point of the plasma.
- calc4DWeights.jl. Compute four-dimensional phase-space weight functions for an (E,p,R,z) grid.
- calcEpRzTopoMap.jl. Compute orbit-type topological maps for a (E,p,R,z) grid. Including maps of poloidal and toroidal transit times, if requested.  
- calcOrbWeights.jl. Compute orbit weight functions, and save them as an orbit weight matrix. Use the DRESS code or simply compute and bin projected velocities if requested.  
- calcSpec.jl. Compute the expected synthetic diagnostic signal spectrum, given a fast-ion and bulk (thermal) plasma distribution.  
- calcTopoMap.jl. Compute orbit-type topological map for a (E,pm,Rm) grid. Include maps of poloidal and toroidal transit times, if requested.  
- ps2os.jl. Transform a fast-ion distribution from (E,p,R,z) space to (E,pm,Rm) space. Save as a 1D compressed vector, and also as a 3D data array if requested.  
- ps2WF.jl. Compute WF signals with abtrirarily high grid resolution in (E,pm,Rm) space. Also compute orbit splits of W, F and WF if requested.
- solveInverseProblem.jl. Given the relationship S=WF, diagnostic measurements S=\[S1,S2,...,Sn\] and corresponding forward model weight matrices W=\[W1,W2,...,Wn\], use tomographic inversion techniques to solve the ill-posed inverse problem and compute the most likely fast-ion distribution F. Collision physics and electromagnetic wave heating in the ion cyclotron range of frequencies (ICRF) physics are supported as prior information.

## 9. OWCF apps
The interactive web applications of the OWCF can be run by navigating to the OWCF/apps/ folder using a command-line terminal and then typing the following:

> `jupyter notebook XWebApp.ipynb`

where X is replaced by the app of your choice.

### 9.1 The apps of the OWCF include:

- comWebApp.ipynb. A simple app that lets the user visualize a Solov'ev magnetic equilibrium and fast-ion orbits. The Solov'ev equilibrium can be customized interactively and the orbit(s) can be changed by changing (E,mu,Pphi; sigma) coordinate values. The app will then re-compute new orbits and magnetic equilibria in real-time.  
- distrWebApp.ipynb. An app that lets the user visualze a (E,pm,Rm) fast-ion distribution in terms of fast-ion energy slices. The user can also choose to input a second fast-ion distribution, and the two will then be able to be compared interactively.  
- EpRzWebApp.ipynb. An app that lets the user visualize orbit topological maps in (E,p,R,z) coordinate space. The option to also visualize poloidal and toroidal transit times exists. The user can also choose to input a fast-ion distribution in (E,p,R,z) format. The app will then superimpose topological boundaries and the fast-ion distribution can be visualzed via a toggle button.  
- modeAnalysisWebApp.ipynb. An app that lets the user visualize simple MHD mode resonances in (E,pm,Rm) orbit space. The user can interactively explore different energies, mode numbers and frequencies. Only toroidal (n) and poloidal (m) mode numbers are included. (E,mu,Phi;sigma) space can be accessed via a toggle button.  
- orbitsWebApp.ipynb. It is arguably the flagship app of the OWCF. It lets the user visualize fast-ion orbits interactively via a topological map in (E,pm,Rm) space. Maps for poloidal and toroidal transit times can also be included. The user can switch to (E,mu,Pphi;sigma) space via a toggle button.  
- orbitWebApp.ipynb. A simple app that lets the user visualize a single fast-ion orbit in detail. The endpoint of the orbit can be changed interactively. The coordinate space is (E,pm,Rm).  
- signalWebApp.ipynb. An app that lets the user visualze a WF signal and orbit splits of WF, W and F, together with their dependence on E, pm and Rm, respectively. Log-scales, splitting, fractions and more options can be changed interactively via toggle buttons.  
- weightsWebApp.ipynb. An app that lets the user interactively visualize the orbit weight functions of an orbit weight matrix. The diagnostic measurement bin and fast-ion energy slice of interest are changed via sliders. A WF signal, an S signal, a fast-ion distribution and null orbits can all be optionally included. The weight visualization can be switched to (E,mu,Pphi;sigma) via a toggle button. Diagnostic viewing cones can be optionally included.  

## 10. Helper scripts
### The helper scripts of the OWCF include:
- calcOrbGrid.jl. Compute an orbit grid and save it, and the valid orbits, as an easily loadable .jld2 file, for the other OWCF scripts to use.  
- calcOrbSpec.jl. Keep track of the functions that compute expected synthetic diagnostic spectra given a fast-ion orbit.
- extractNullOrbits.jl. Take a diagnostic signal and an orbit weight matrix in 2D form, and identify the null orbits, given thresholds.  
- extractTempDens.jl. Take TRANSP data and extract the thermal plasma density and temperature profiles.  
- extractTopoBounds.jl. Take a topological map and examine it to locate the boundaries between topological regions. Save the topological boundaries as a .jld2 file.  
- F_os_1Dto3D.jl. Inflate a fast-ion distribution in its compressed 1D form into its full 3D form.  
- orbweights2Dto4D.jl. Inflate an orbit weight matrix from its compressed 2D form into its full 4D form, set by the user.  
- orbWeights4Dto2D.jl. Take an inflated orbit weight matrix in its 4D form, and compress it into its 2D form.  
- os2com.jl. A script that can be used to transform quantities from (E,pm,Rm) orbit space to (E,mu,Pphi;sigma) constants-of-motion space.  

## 11. Extra scripts
### The extra scripts of the OWCF include:  
- constants.jl - A file that keeps track of all the constants used by various parts of the OWCF
- createCustomFIDistr.jl. Lets the user create a custom fast-ion distribution function to be used as input in other OWCF scripts.
- createCustomLOS.jl. This script will help the user create a custom line-of-sight (LOS) file, usable by the OWCF. The file format mimics the file format used by the LINE21 code. The options for the custom LOS includes detector position, angle and size. A cylinder shape is assumed for the LOS.
- createCustomMagneticEquilibrium.jl. Compute a Solov'ev magnetic equilibrium and first wall geometry, based on the inputs. Save as .jld2 file.
- createCustomThermalProfiles.jl. Lets the user create custom thermal plasma species profiles (temperature and density), to be used as input in other OWCF scripts.
- dependencies.jl. A collection of many functions used by various parts of the OWCF.
- gui.jl. A collection of functions that can be used to plot various quantities computed by the OWCF.
- TRANSP_fastion_data_to_EpRz.jl. Given a TRANSP fast-ion distribution file and user-requested (E,p,R,z) grid, compute the fast-ion distribution on that grid. Save as .jld2 file.

## 12. Misc scripts
### The misc scripts of the OWCF include:
- availReacts.jl. Keep track of the fusion reactions availabe for simulation by the OWCF, through the DRESS code.
- convert_units.jl. A collection of functions to work with units of measurements, e.g. compute conversion factors between different units of measurements.
- diag2tokamak.jl. A (collection of) function(s) that keeps track of which diagnostics are installed at which tokamaks.
- load_TRANSP_interp_object. A script that contains a function to be able to load TRANSP data as interpolation objects. These objects are then used in various parts of the OWCF.
- rewriteReacts.jl. A collection of functions that provides the utility to re-write the notation of fusion reactions in different formats.  
- species_func.jl. A collection of functions that provide easy access to the properties of various ion species.  
- temp_n_dens.jl. Store and provide default thermal plasma temperature and density profiles for the OWCF. Also includes loading and interpolation functions for TRANSP thermal temperature and density data.  

Other data in the misc/ folder includes: default_temp_n_dens.png, eqdsk_file_breakdown.pdf, howToLoadFromTRANSP_FI_CDF.jl, logo.gif and transp_outputs.txt.

## 13. OWCF on HPC clusters
> If you would like to use the OWCF on a SLURM computational cluster, you can use the .sh submit file templates to easily submit OWCF runs for batch jobs. Templates for other computational cluster workload managers than SLURM are currently not included with the OWCF.

## Good luck! Please consult the OWCF_manual.pdf, howToInstallJuliaAndTheOWCF_MacOS.pdf and howToInstallJuliaAndTheOWCF_Windows.pdf documents for further info.

Henrik Järleblad
A (tired again) postdoc
October 9th, 2025