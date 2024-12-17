# This repository contains useful tools for the molecular dynamics package LAMMPS

## Python-for-LAMMPS:
Scripts to handle LAMMPS output files.

## coulomb_smooth_potential:
Scripts to perfom an additional type of potential that allows the computation of the Coulomb interaction using a smooth truncation in LAMMPS.

## TersoffModZbl:
Scripts to perfom an additional type of potential that allows the computation of the Modified Tersoff potential joined to the ZBL potential in LAMMPS.

## TempVcut:
Scripts to compute the temperature estimator avoiding atoms that moves above a threshold velocity
This compute calculates the coarsed-grained ionic temperature estimator in the two-temperature model according to the cutoff approach (https://dx.doi.org/10.1088/1361-648X/ad4941).

## FixTTM_MMMG:
Scripts to introduce a reformulation of the traditional 2T-MD scheme in order to overcome the limitations on the treatment of the high-energy particles by addressing the spurious double-interaction of fast-moving atoms with electrons.
This fix switches off the electron-phonon terms in fast-moving particles that exceed a threshold velocity, only applying electronic stopping effects.
This corresponds to the cutoff approach, explained in (https://dx.doi.org/10.1088/1361-648X/ad4941).