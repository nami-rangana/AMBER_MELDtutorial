# AMBER INSTALLATION GUIDE
This guide provides installation instructions for AMBER on the HiPerGator research computing cluster as a locally installed package.

## 1. open session:
Execute the following command on a HiPerGator login node to access a compute node. This requests a single node with 1 GPU for a 5-hour period:

```srun --nodes=1 --ntasks=4 --gpus=1 -m 32 -t 5:00:00 --pty /bin/bash -i```

## 2. load required modules:
Now we need to load the following modules:
<!-- lapack/3.11.0 cmake/3.30.5 -->
```ml cuda/12.8.1 gcc/14.2.0 openmpi/5.0.7```

## 3. Extract downloaded setup files:
Navigate to your installation directory and extract the source archives.

#### AmberTools: 
```
cd <amber_path>
tar xvjf ambertools25.tar.bz2
```
#### Amber (pmemd): 
```
cd <amber_path>
tar xvjf pmemd24.tar.bz2
```

## 4. Replace edited files -- MELD _(optional)_
If you have modified source code, update the files before building:

```
rm -f <amber_path>/pmemd24_src/src/pmemd/src/<files>
cp <my_files> <amber_path>/pmemd24_src/src/pmemd/src/.
```

## 5. Build and Install
#### AmberTools: 
If you have tried installing before you might still have some old files. Just to be safe we can remove them. (optional)
```
cd <amber_path>/ambertools25_src/
rm -rf build && mkdir build && cd build
```
Then copy [run_cmake](run_cmake) file to the build directory and run it with `chmod +x run_cmake; ./run_cmake`. Alternatively use following command:
```
cmake <amber_path>/ambertools25_src \
    -DCMAKE_INSTALL_PREFIX=<amber_path>/ambertools25 \
    -DCOMPILER=GNU  \
    -DMPI=TRUE -DCUDA=TRUE -DINSTALL_TESTS=TRUE \
    -DDOWNLOAD_MINICONDA=TRUE \
    2>&1 | tee  cmake.log
```
**Note:** Specify your desired installation directory for `-DCMAKE_INSTALL_PREFIX=`.<br>If errors are reported, search for 'CMake Error' in the cmake.log file.
If no errors were encountered, compile and install:
```
make install -j4
```
**Note:** Use `-j1` if errors occur during compilation.

Now this will install AmberTools25 to `<amber_path>/ambertools25`. The installation will create a resource file `amber.sh`. Source this file to configure your shell environment for AMBER:
```
source <amber_path>/ambertools25/amber.sh
```

#### Amber (pmemd): 
If you have attempted installation previously, remove old build files to ensure a clean build (optional):
```
cd <amber_path>/pmemd24_src/
rm -rf build && mkdir build && cd build
```
Configure the build using CMake:
```
cmake $AMBER_PREFIX/pmemd24_src \
    -DCMAKE_INSTALL_PREFIX=$AMBER_PREFIX/pmemd24 \
    -DCOMPILER=GNU  \
    -DMPI=TRUE -DCUDA=TRUE -DINSTALL_TESTS=TRUE \
    -DDOWNLOAD_MINICONDA=TRUE -DBUILD_PYTHON=FALSE \
    -DBUILD_PERL=FALSE -DBUILD_GUI=FALSE \
    -DPMEMD_ONLY=TRUE -DCHECK_UPDATES=FALSE \
    2>&1 | tee  cmake.log
```
**Note:** Specify your desired installation directory for `-DCMAKE_INSTALL_PREFIX=`.<br>If errors are reported, search for 'CMake Error' in the cmake.log file.
If no errors were encountered, compile and install:
```
make install -j4
```
**Note:** Use `-j1` if errors occur during compilation.

Now this will install Amber to `<amber_path>/pmemd24`. The installation will create a resourcee file `amber.sh`. This file will set up your shell environment correctly for Amber after sourcing.
```
source <amber_path>/pmemd24/amber.sh
```
<!--
## 5. configure for L4 (Ada, sm_89) and reduce ptxas optimization using FindCUDA vars

```
cmake -S .. -B . \
  -DCOMPILER=GNU \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/orange/alberto.perezant/AMBER-MELD/amber25 \
  -DMPI=ON -DCUDA=ON -DOPENMP=ON -DDOWNLOAD_MINICONDA=ON -DINSTALL_TESTS=ON \
  -DCUDA_ARCH_NAME=Manual -DCUDA_ARCH_BIN=89 -DCUDA_ARCH_PTX=89 \
  -DCUDA_NVCC_FLAGS_RELEASE="-Xptxas=-O0 --maxrregcount=48"
```
##  6. build the CUDA target strictly serial
This step makes it easy to troubleshoot

```cmake --build . --target pmemd_cuda_SPFP -j1```

## 7. build the rest
```cmake --build . -j6```

## 8. install
```cmake --install .```
-->
