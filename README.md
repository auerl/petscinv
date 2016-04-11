PETScinv
========

```fortran
                                     _|                          _|
               _|_|_|      _|_|    _|_|_|_|    _|_|_|    _|_|_|      _|_|_|    _|      _|
               _|    _|  _|_|_|_|    _|      _|_|      _|        _|  _|    _|  _|      _|
               _|    _|  _|          _|          _|_|  _|        _|  _|    _|    _|  _|
               _|_|_|      _|_|_|      _|_|  _|_|_|      _|_|_|  _|  _|    _|      _|
               _|
               _|
```

PETScinv solves global tomographic imaging problems in parallel using the PETSc component KSP. The scalable linear equations solvers (KSP) module provides an easy-to-use interface to the combination of a Krylov subspace iterative method and a preconditioner (in the KSP and PC components, respectively) or a sequential direct solver. 

The code currently only supports orthogonal, curvilinear hexahedral basis functions (spherical sections or voxels), but adding modules to treat other parametization strategies (splines, spherical harmonics, wavelets) to discretize the model space, should be easy to add. One can account for different regularization strategies, such as roughness minimization or norm minimization. Tomography results can be exported in a simple ASCII format or as netcdf files, together with an .xdmf file, ready to be explored with [ParaView](http://www.paraview.org/).
  
Quick start
===========

The code is easy to run. Install on (fat memory) compute node using

```bash
git clone http://github.com/auerl/petscinv.git
```

Make sure that PETSc library is installed (on HPC systems this is often a standard library that you can load via the module environment and commands like `module load petsc`. We have tested the code with GCC and the Intel Compiler Collection. Compile the code with `make`. Alongside with the code comes some ray-theoretical tomography submatrices for a smallish dataset. Inverting this dataset using a parallel implementation of the GMRES algorithm can by done by typing `run_petscinv`, or commands like

```bash
mpirun -n 100 ./petscinv \
    -matrix_schedule 'samples/schedule' \
    -inparam_file 'samples/inparam' \
    -inversion_parameters 'vsh,vsv' \
    -reference_model 'prem_ani' \
    -number_of_layers 21 \
    -equatorial_increment 5.0 \
    -project_id 'INV01' \
    -output_format 'ascii' \
    -solution_type 'normal' \
    -ksp_monitor \
    -grouped_varr \
    -ksp_type gmres
```

Please see the documentation for our toolbox [flexinv](../flexinv) to find out how to assemble ray-theoretical or waveform based tomography matrices for realistic imaging problems.
