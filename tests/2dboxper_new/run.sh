cd ../../
make uninstall
make install

cd tests/2dboxper_new/
rm -rf vtk/
rm -rf obj/
./setup maxwell 2dboxper_new
make
#./nekcem
