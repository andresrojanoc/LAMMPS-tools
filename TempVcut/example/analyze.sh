#!/bin/bash

#Script to analyze the results

#Extract the Ionic grids with standard temperature
python3 grid.py CIT.out > temp_CIT.out

#Extract the Ionic grids with cutoff temperature
python3 grid.py CIT_cut.out > temp_CIT_cut.out

#Delte temporary grids
rm *grid_.out

#Delte CIT grids (comment to keep them)
rm CIT_*0.out