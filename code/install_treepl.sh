#!/bin/bash

# Install treePL
# Set some important environment variables
export BOOST_ROOT=$CONDA_PREFIX
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${CONDA_PREFIX}/lib64:${LD_LIBRARY_PATH}"
export C_INCLUDE_PATH=$CONDA_PREFIX/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$CONDA_PREFIX/include:$CPLUS_INCLUDE_PATH
export DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib:$DYLD_LIBRARY_PATH

# Download and uncompress
wget https://github.com/blackrim/treePL/archive/refs/heads/master.zip
unzip master.zip && rm master.zip && mv treePL-master treePL

# Install dependencies, beginning with nlopt
cd treePL/deps
tar -xzvf nlopt-2.4.2.tar.gz
cd nlopt-2.4.2
# Make sure the configuration script points to the conda environment
./configure --prefix=$CONDA_PREFIX LDFLAGS="-L${CONDA_PREFIX}/lib -fPIC" CFLAGS="-I${CONDA_PREFIX}/include -fPIC"
make
make install

# And then ADOL-C This is a bit more sensitive.
cd ../
tar -xzvf ADOL-C-2.6.3.tgz
cd ADOL-C-2.6.3
# Make sure the configuration script points to the conda environment and conda-installed boost libraries
./configure --prefix=$CONDA_PREFIX --with-openmp-flag=-fopenmp LDFLAGS="-L${CONDA_PREFIX}/lib64" CFLAGS="-I${CONDA_PREFIX}/include" --with-boost=$CONDA_PREFIX --with-boost-libdir=$CONDA_PREFIX/lib --with-boost-include=$CONDA_PREFIX/include
make
make install

# Finally, install treePL
cd ../../src
# And configure. NOTE: This configure script does not recognize provided LDFLAGS or CFLAGS. It hard codes them to /usr/local. We will need to modify the Makefile directly.
./configure --prefix=$CONDA_PREFIX
# Update them to work with conda:
sed "s|/usr/local|${CONDA_PREFIX}|g" Makefile > tmp && mv tmp Makefile
sed "s|/usr/lib64:||g" Makefile > tmp && mv tmp Makefile
sed "s|/usr/bin/|${CONDA_PREFIX}/bin/|g" Makefile > tmp && mv tmp Makefile
make
make install
# remove the temporary directory used for building treePL and PATHd8
cd ../../ && rm -r tmpdir
