# AxiSEM 1.4 [![Build Status](https://travis-ci.org/geodynamics/axisem.svg?branch=master)](https://travis-ci.org/geodynamics/axisem)

## Axially symmetric Spectral Element Method

Copyright 2018, Tarje Nissen-Meyer, Martin van Driel, Simon Stähler, Kasra Hosseini, Stefanie Hempel, Lion Krischer, Alexandre Fournier

Webpage and distribution: http://www.axisem.info
Contact and information:  info@axisem.info

April, 13, 2018 

## Citation
If you are publishing results obtained with this code, please cite this paper:

T. Nissen-Meyer, M. van Driel, S. C. Staehler, K. Hosseini, S. Hempel, L. Auer, A. Colombi and A. Fournier:
**"AxiSEM: broadband 3-D seismic wavefields in axisymmetric media"**, *Solid Earth*, 5, 425-445, 2014
doi:10.5194/se-5-425-2014 http://www.solid-earth.net/5/425/2014/

## Content of the repository
`manual_axisem_1.4.pdf` - PDF manual

`MESHER` - The program to generate 2D meshes for the SEM forward solver

`SOLVER` - the SEM forward solver itself

`submit.py` - a Python script to create an Instaseis database

`make_axisem.macros` - macro file to set compiler options

`copytemplates.sh` - reset all input files to default templates 

`COPYING` - The GNU General Public License

`HISTORY` - changelog

`README` - this file

## Basic instructions for running:

More details are found in the manual. For a quick start:

 - The mesher has to be run first to generate a SEM mesh for the solver. 
 - Any changes in the resolution, spherically symmetric background model, or number 
   of processors requires a new mesh. 
 - Changes in the source-receiver settings or 2.5D heterogeneities only need a new solver run.
 - Changes on the moment tensor, receiver components, or filtering only need a new postprocessing run.
 - Using Instaseis, seismograms for any depth or moment tensor can be calculated from the wavefield of one force source at the surface.

General settings and explanations for parameters needed in MESHER and SOLVER 
are found in the `inparam_*` files in the respective directories. 

To create Instaseis databases quickly, run the script submit.py in the main directory.

1) Run `copytemplates.sh` to set up a generic run with pre-set parameters

2) Go into the MESHER directory, run `./submit.csh`. This compiles the code using
gfortran and mpif90 as default compilers, and then submits a job on a single node. 
For high resolution (seismic period below 3s), this requires significant amounts 
of RAM, see manual.

3) Check `OUTPUT`; if finished, then run `./movemesh.csh <MESH_NAME>`

4) go into `SOLVER`, check the vtk files in `/MESHES/<MESH_NAME>` with Paraview, if you want

5) Edit `inparam_basic` to the desired `<MESH_NAME>`

6) Run `./submit.csh <RUN_NAME>` , this compiles and then submits a parallel run.

7) Check `<RUN_NAME>/OUTPUT*`, which will record the progress of the run and any problems.

8) All data-related output is in `<RUN_NAME>/Data/`

For `STATIONS` receivers, the solver places each station at depth
`burial_depth - elevation` (metres), clipped between zero and the model's
outer radius. It searches the elastic mesh at that depth and interpolates
displacement within the containing element. A station in the fluid domain or
outside the mesh cannot be recorded.

### Boundary-face wavefields

Set `KERNEL_WAVEFIELDS true`, `USE_NETCDF true`, and `SAVE_BDRY_FACES true`
in `SOLVER/inparam_advanced`. Generate `SOLVER/boundary_faces.dat` with
AxiSEMLib's `surface_merge.py`, or edit `SOLVER/boundary_surfaces.dat` to contain your
face coordinates. `submit.csh` copies the selected input into the run directory. The file
ships with five example faces: four solid (100, 500, 1500, and 2500 km) and
one fluid (4000 km). The file has no header and contains one
`longitude latitude depth_km` row per point,
ordered as consecutive 5×5 faces. Each face's 13th row determines whether its
25 points are sampled in the solid or fluid; points outside that phase cause
an error. `DUMP_T0` and `KERNEL_SPP` control the sampling times.

The run writes `Data/boundary_wavefields.nc4`. It stores each distinct element
once: `disp_s` and `disp_z` have dimensions `(time, solid_element, npol, npol)`,
`disp_p` has the same dimensions for dipole and quadrupole sources, and `chi`
has dimensions `(time, fluid_element, npol, npol)`. The last two axes are the
element's GLL nodes. The displacement components retain AxiSEM's source-centered
cylindrical modal convention.

`boundary_wavefields.nc4/Mesh` contains geometry, material properties, GLL
coordinates, and connectivity for these distinct elements. The 0-based
`solid_element_to_mesh` and `fluid_element_to_mesh` arrays index `Mesh/elements`. The
`point` dimension retains the input faces for lookup: `face_phase` is 0 for
solid or 1 for fluid, `point_element_slot` selects the appropriate wavefield
element, and `xi`/`eta` locate the original point within it. `element_id` is
the 0-based global solver ID for each point. `source_phi` gives its azimuth in
the source frame. The regular kernel wavefields and Mesh group are omitted
from `axisem_output.nc4` in this mode.
The solver log reports each point's owning process, global element, and
meridional location error (m).

9) Convolution with a STF and summation for a moment source is done by running
   `postprocessing.csh` in `<RUN_NAME>`

10) A more modern and efficient way of seismogram retrieval is using Instaseis,
    a Python toolbox to retrieve seismograms for arbitrary depths and moment
    tensors from the stored wavefield of one AxiSEM run (http://www.instaseis.net)

Detailed instructions can be found in the file `manual_axisem_1.4.pdf`

## Installation
Generally, AxiSEM needs
- a C and a Fortran90 compiler (typically gcc and gfortran)
- MPI (typically OpenMPI)
- NetCDF (not absolutely necessary, but needed to create an Instaseis database)

### MacOS X
The easiest way to install the necessary requirements is *homebrew* (https://brew.sh/). If you have not installed homebrew itself yet, do so as described on the homepage. The preferred way (as of May 2022) is 
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
and follow the steps therein.

Afterwards, install gfortran, openmpi and NetCDF with 
```
brew install gfortran gcc openmpi netcdf hdf5-mpi
```

### Linux
The necessary packages can easily be installed using APT
```
sudo apt install gfortran gcc libopenmpi-dev openmpi-bin libnetcdff-dev 
```
