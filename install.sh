# compile divsufsort
mkdir -p libdivsufsort-master/build
cd libdivsufsort-master/build
cmake -DCMAKE_BUILD_TYPE="Release" -DBUILD_DIVSUFSORT64="ON" -DBUILD_SHARED_LIBS="OFF" -DCMAKE_INSTALL_PREFIX=`pwd` ..
make
make install
cd ../../

# make GeDi
make
