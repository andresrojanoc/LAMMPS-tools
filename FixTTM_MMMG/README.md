# Script to modify the two-temperature model (ttm) in LAMMPS

These changes are to introduce a reformulation of the traditional 2T-MD scheme to overcome the limitations on the treatment of the high-energy particles by addressing the spurious double-interaction of fast-moving atoms and electrons.

(https://dx.doi.org/10.1088/1361-648X/ad4941)

## Modifications:

(1) Ensure that high-energy particles are subject to electronic stopping only and not to electron-phonon coupling interactions (as this is an equilibrium process).

```
keyword: cutoff
```

(2) Electronic specific heat is given using a table file.

```
keyword: table <file name>
```

## Content:

**fix:** Directory with the scripts with the reformulated ttm, one with the original fix FixTTM (ttm) and an additional FixTTM_MMMG (ttm_mmmg).

**example:** Directory with an example of the reformulated ttm.

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

Include the desired LAMMPS' modules (e.g., MANYBODY, EXTRA-FIX)

```
make yes-manybody
```

```
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

To extract the ionic temperature appropriately, use TempVcut compute. (check the correct ./path/)

```
cp ./path/compute_temp_vcut_chunk.* .
```

Compile LAMMPS

```
make mpi #(paralel)
```

```
make serial #(serial)
```

## Notes:

See the in.pka file in the example directory for implementation details and appropriate usage.

Extract the configuration file in the example directory using tar -xvzf equi.tar.gz.

The table file interpolates the values of the electronic specific heat for each electronic cell.

The values in the table are provided in terms of (temperature, first column) and (energy/volume/temperature, second column). For example: units of K and eV / (A^3 · K).

The CIT_cut compute does not print real temperatures when the cutoff key is disabled (default).

The files grid.py, gridE.py, and analyze.sh are provided as tools to assist with post-processing.