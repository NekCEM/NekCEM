cd ../../
make uninstall
make install

cd tests/3dboxper_new/
rm -rf vtk/
rm -rf obj/
rm nekcem
./setup maxwell 3dboxper_new
make
#./nekcem
