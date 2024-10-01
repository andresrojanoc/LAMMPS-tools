## Script to perfom an additional type of compute in LAMMPS
## This new compute calculates the temperature neglecting (athermal) fast-moving atoms

compute: Directory with the scripts to compile ComputeTempChunkVcut (temp_vcut/chunk)
example: Directory with an example of an input file that employs this compute

## Ussage:

# Extract the downloaded tarball (https://www.lammps.org/download.html) in the desired directory

tar -xvzf lammps-stable.tar.gz

# Load the module (if necessary)

module load iomkl

# Change directory to the STUBS directory in the src

cd lammps/src/STUBS

# Perform the make command

make

# Go back to the src

cd ../
(check with pwd)

# Clean all possible residues (perhaps this is not necessary when a new LAMMPS is compiled)

make clean-all

# Include the desired LAMMPS' modules (e.g., MANYBODY, EXTRA-FIX)

make yes-manybody
make yes-extra-fix

# Copy the content of the compute folder (i.e., compute_temp_vcut_chunk.h and compute_temp_vcut_chunk.cpp) to the /src/

cp ../../compute/* .
(check the correct path where you saved the smooth directory)

#Compile LAMMPS

make mpi (paralel)
make serial (serial)

Notes:

See in.pka in the example directory for the implementation and appropriate usage
The files compute_temp_vcut_chunk are based on compute_temp_chunk
The example is just for showing purposes, as the ttm is the standard and not the cutoff approach
In other words, the example is just to showcase the functionality of the compute, neglet fast-moving atoms in the ionic temperature estimator
