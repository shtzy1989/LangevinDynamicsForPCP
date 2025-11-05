rm -rf build
cmake -B build -DCMAKE_INSTALL_PREFIX=install
cd build
make -j 12
make install
cd ../
