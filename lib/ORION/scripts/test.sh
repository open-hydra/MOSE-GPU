#!/bin/bash

cd $ORIONDIR/bin/test
rm -f XML* *.dat *.plt field*

ulimit -s unlimited

#echo
#echo '--- ORION ---'
#./vtk_fullpower -strgr
#ORION --out-format=ascii --in-format=raw

echo
echo '--- VTK writing (wrapper) ---'
./vtk_write_wrapper

echo
echo '--- VTK reading (wrapper) ---'
./vtk_read_wrapper

echo
echo '--- Tecplot writing ---'
./tecplot_write

echo
echo '--- Tecplot reading (ascii) ---'
./tecplot_read


echo
echo '--- PLOT3D writing ---'
./p3d_write

echo
echo '--- PLOT3D reading ---'
./p3d_read
