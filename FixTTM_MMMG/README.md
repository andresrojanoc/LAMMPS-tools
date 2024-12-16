# Script to modify the two-temperature model (ttm) in LAMMPS

These changes are to introduce a reformulation of the traditional 2T-MD scheme to overcome the limitations on the treatment of the high-energy particles by addressing the spurious double-interaction of high-energy atoms with electrons.

(https://dx.doi.org/10.1088/1361-648X/ad4941)

## Modifications:

(1) Removed the thermalisation time from the start of the cascade by ensuring that that electronic stopping is only applied
to the high-energy particles and not electron-phonon coupling (as this is an equilibrium process)

(2) If a particle is subject to electronic stopping it is not included in the temperature for the ionic voxel, to avoid
fast-moving particles causing the temperature of a whole voxel to be very high leading to electron-phonon coupling that should not occur

(3) Additionally, the electronic specific heat is given the option to be called through a Table


## Content:

fix: Directory with the scripts with the reformulated ttm, with the original fix FixTTM (ttm) or an additional FixTTM_MMMG (ttm_mmmg)

example: Directory with an example of the reformulated ttm

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

Go back to the src
```
cd ../
```
(check with pwd)

Clean all possible residues (perhaps this is not necessary when a new LAMMPS is compiled)
```
make clean-all
```

Include the desired LAMMPS' modules (e.g., MANYBODY, EXTRA-FIX)

```
make yes-manybody
make yes-extra-fix
```

You can either replace the original fix by replacing its files (i.e., fix_ttm.h and fix_ttm.cpp) in the folder /src/EXTRA-FIX/ (check the correct path)

```
cp ../../fix/fix_ttm.* EXTRA-FIX/.
```

Or you can keep the original and create an additional fix (i.e., fix_ttm_mmmg.h and fix_ttm_mmmg.cpp) in the folder /src/
```
cp ../../fix/fix_ttm_mmmg.* .
```

To extract the ionic temperature appropriately, use TempVcut compute.
```
cp ./path/compute_temp_vcut_chunk.* .
(change ./path/, with the folder where the files for the ComputeTempChunkVcut is located)
```

Compile LAMMPS

```
make mpi (paralel)
make serial (serial)
```

## Notes:

See the in.pka file in the example directory for implementation details and appropriate usage

Extract the configuration file in the example directory using tar -xvzf equi.tar.gz

The table file interpolates the values of the electronic specific heat for each electronic cell

The values in the table are provided in terms of (temperature, first column) and (energy/ volume * temperature, second column). For example: units of K and eV / (AA^3 · K)

The CIT_cut compute does not print real temperatures when the cutoff key is disabled (default behavior)

The files grid.py, gridE.py, and analyze.sh are provided as tools to assist with post-processing