# Script to perfom an additional type of Tersoff potential in LAMMPS

This new potential allows to compute the Modified Tersoff joined to the ZBL.

## Content:

**pair:** Directory with the scripts to compile LAMMPS with the TersoffModZbl interaction (tersoff/mod/zbl).

**example:** Directory with an example of an input file that employs this potential.

## Ussage:

Extract the downloaded tarball (https://www.lammps.org/download.html) in the desired directory

```
tar -xvzf lammps-stable.tar.gz
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

Include the desired LAMMPS' modules (e.g., MANYBODY)

```
make yes-manybody
```

Copy the content of the pair folder (i.e., pair_tersoff_mod_zbl.cpp, pair_tersoff_mod_zbl.h) to the /src/ (check the correct path where you saved the smooth directory)

```
cp ../../pair/* .
```

Compile LAMMPS

```
make mpi #(paralel)
```

```
make serial #(serial)
```
Notes

See in.file in the example directory for the implementation.

The first 20 terms of the potential file are those explained in the pair_tersoff_mod. The last 4 terms of the potential file are the last four of the pair_tersoff_zbl.
