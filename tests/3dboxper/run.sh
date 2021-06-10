cd ../../
make uninstall
make install

cd tests/3dboxper/
rm -rf vtk/
rm -rf obj/
rm nekcem
./setup maxwell 3dboxper
make
#./nekcem
