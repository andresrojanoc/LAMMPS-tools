# Scripts to perfom an additional type of Coulomb potential in LAMMPS

This potential allows to compute the Coulomb interaction using a smooth truncation.

## Content:

**smooth:** Directory with the scripts to compile LAMMPS with the Coulomb interaction using a smooth truncation (coul/smooth).

**example:** Directory with an example of an input file that employs this potential and the corresponding configuration file.

## Ussage:

Extract the downloaded tarball (https://www.lammps.org/download.html) in the desired directory

```
tar -xvzf lammps-stable.tar.gz.
```

Load the module (if necessary)

```
module load iomkl
```

Change directory to the STUBS directory in the src

```
cd lammps/src/STUBS
```

Perform the make command

```
make
```

Go back to the src (check with pwd)

```
cd ../
```

Clean all possible residues

```
make clean-all
```

Include the desired LAMMPS' modules (e.g., molecule, kspace, rigid)

```
make yes-molecule
```

Copy the content of the smooth folder (i.e., pair_coul_smooth.cpp, pair_coul_smooth.h) to the /src/ (check the correct path where you saved the smooth directory)

```
cp ../../smooth/* .
```

#Compile LAMMPS

```
make mpi #(paralel)
```

```
make serial #(serial)
```

Notes:

The cutoff radius is set on the selected units.
See in.smooth in the example directory for the implementation.
This potential corresponds to that explained in Walther et al. (2001) J. Phys. Chem. B DOI: 10.1021/jp011344u.
