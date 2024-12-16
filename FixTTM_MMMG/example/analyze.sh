#!/bin/bash

#Script to analyze the results

echo "Extract ionic grid temperature with standard approach"
python3 grid.py CIT.out > Ta.temp

echo "Extract ionic grid temperature with cutoff approach"
python3 grid.py CITcut.out > Tacut.temp

echo "Extract electronic grid temperature"
python3 gridE.py CET.out.1 100 100 > Te.temp

echo "Delete temporary grids"
rm *grid_.out

echo "Delete CIT grids"
rm CIT_*.out
rm CITcut_*.out

echo "Move electronic temperatures to a E_grid folder"
mkdir E_grid
mv CET.out* E_grid
