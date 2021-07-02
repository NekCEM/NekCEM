CEED_YIMIN="$HOME/Desktop/github/ceed-yimin/codes/2dboxper"
CASE_NAME="2dboxper"
CASENAME=$CASE_NAME
CASE_NAME_NEW="2dboxper_new"
export CASENAME

cd ../../
sudo make uninstall_inplace
sudo make install

cd tests/2dboxper/

cp $CEED_YIMIN/$CASE_NAME.usr ./ 
sudo rm -rf vtk/
sudo rm -rf obj/
sudo rm -rf logs/*
sudo rm nekcem
#sudo ./setup maxwell 2dboxper
configurenek --FC mpif77 maxwell 2dboxper
./mkuserfile
sudo make

###########################
####### .map + .ma2 #######
###########################
rm -rf $CASE_NAME.*

# old .rea + .map
cp $CEED_YIMIN/$CASE_NAME.rea ./ 
cp $CEED_YIMIN/$CASE_NAME.map ./
./nekcem &> logs/output_rea_map
echo '===== Finish .rea+.map serial run'
mpirun -np 2 ./nekcem &> logs/output_rea_map_mpinp2
echo '===== Finish .rea+.map mpi run'
rm -rf $CASE_NAME.*

# old .rea + .ma2
cp $CEED_YIMIN/$CASE_NAME.rea ./ 
cp $CEED_YIMIN/$CASE_NAME_NEW.ma2 ./$CASE_NAME.ma2
./nekcem &> logs/output_rea_ma2
echo '===== Finish .rea+.ma2 serial run'
mpirun -np 2 ./nekcem &> logs/output_rea_ma2_mpinp2
echo '===== Finish .rea+.ma2 mpi run'
rm -rf $CASE_NAME.*

# new .rea + .re2 + .map
cp $CEED_YIMIN/$CASE_NAME_NEW.rea ./$CASE_NAME.rea
cp $CEED_YIMIN/$CASE_NAME_NEW.re2 ./$CASE_NAME.re2
cp $CEED_YIMIN/$CASE_NAME.map ./$CASE_NAME.map
./nekcem &> logs/output_rea_re2_map
echo '===== Finish .rea+.re2+.map serial run'
mpirun -np 2 ./nekcem &> logs/output_rea_re2_map_mpinp2
echo '===== Finish .rea+.re2+.map mpi run'
rm -rf $CASE_NAME.*

# new .rea + .re2 + .ma2
cp $CEED_YIMIN/$CASE_NAME_NEW.rea ./$CASE_NAME.rea
cp $CEED_YIMIN/$CASE_NAME_NEW.re2 ./$CASE_NAME.re2
cp $CEED_YIMIN/$CASE_NAME_NEW.ma2 ./$CASE_NAME.ma2
./nekcem &> logs/output_rea_re2_ma2
echo '===== Finish .rea+.re2+.ma2 serial run'
mpirun -np 2 ./nekcem &> logs/output_rea_re2_ma2_mpinp2
echo '===== Finish .rea+.re2+.ma2 mpi run'
rm -rf $CASE_NAME.*

# .par + .re2 + .map
cp $CEED_YIMIN/$CASE_NAME_NEW.par ./$CASE_NAME.par
cp $CEED_YIMIN/$CASE_NAME_NEW.re2 ./$CASE_NAME.re2
cp $CEED_YIMIN/$CASE_NAME.map ./$CASE_NAME.map
./nekcem &> logs/output_par_re2_map
echo '===== Finish .par+.re2+.map serial run'
mpirun -np 2 ./nekcem &> logs/output_par_re2_map_mpinp2
echo '===== Finish .par+.re2+.map mpi run'
rm -rf $CASE_NAME.*

# .par + .re2 + .ma2
cp $CEED_YIMIN/$CASE_NAME_NEW.par ./$CASE_NAME.par
cp $CEED_YIMIN/$CASE_NAME_NEW.re2 ./$CASE_NAME.re2
cp $CEED_YIMIN/$CASE_NAME_NEW.ma2 ./$CASE_NAME.ma2
./nekcem &> logs/output_par_re2_ma2
echo '===== Finish .par+.re2+.ma2 serial run'
mpirun -np 2 ./nekcem &> logs/output_par_re2_ma2_mpinp2
echo '===== Finish .par+.re2+.ma2 mpi run'
rm -rf $CASE_NAME.*

#####################
####### .co2  #######
#####################
cp $CEED_YIMIN/$CASE_NAME.usr ./ 
sudo rm -rf vtk/
sudo rm -rf obj/
configurenek --FC mpif77 --extra-FFLAGS ' -DPARRSB -DPARMETIS' --extra-CFLAGS '-DPARRSB -DPARMETIS' maxwell 2dboxper
./mkuserfile
sudo make

# new .rea + .re2 + .co2
cp $CEED_YIMIN/$CASE_NAME_NEW.rea ./$CASE_NAME.rea
cp $CEED_YIMIN/$CASE_NAME_NEW.re2 ./$CASE_NAME.re2
cp $CEED_YIMIN/$CASE_NAME_NEW.co2 ./$CASE_NAME.co2
./nekcem &> logs/output_rea_re2_co2
echo '===== Finish .rea+.re2+.co2 serial run'
mpirun -np 2 ./nekcem &> logs/output_rea_re2_co2_mpinp2
echo '===== Finish .rea+.re2+.co2 mpi run'
rm -rf $CASE_NAME.*

# .par + .re2 + .co2
cp $CEED_YIMIN/$CASE_NAME_NEW.par ./$CASE_NAME.par
cp $CEED_YIMIN/$CASE_NAME_NEW.re2 ./$CASE_NAME.re2
cp $CEED_YIMIN/$CASE_NAME_NEW.co2 ./$CASE_NAME.co2
./nekcem &> logs/output_par_re2_co2
echo '===== Finish .par+.re2+.co2 serial run'
mpirun -np 2 ./nekcem &> logs/output_par_re2_co2_mpinp2
echo '===== Finish .par+.re2+.co2 mpi run'
rm -rf $CASE_NAME.*



# View output
tail logs/*
