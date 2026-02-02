# AMBER INSTALLATION GUIDE
This file contains installation steps for HiperGator research computing cluster as a local instaled package.
## 1. open session:
Run following command on a Hipergator login node to get access to a compute node. Here we request a single node with 1 GPU for 5 hr time period.

```srun --nodes=1 --ntasks=14 --gpus=1 -m 32 -t 5:00:00 --pty /bin/bash -i```

## 2. load required modules:
Now we need to load the following modules.

```ml cuda/12.8.1 gcc/14.2.0 openmpi/5.0.7 lapack/3.11.0 cmake/3.30.5```

## 3. remove the build directory and create new _(optional)_
If you have tried installing before you might still have some old files. Just to be safe we can remove them.

```
cd <amber_path> #/orange/alberto.perezant/AMBER-MELD/amber
rm -rf build && mkdir build && cd build
```
## 4. replace the edited files _(optional)_
Incase you have edited source code, make sure to update the files.

```
rm -f <files>
cp <my_files> .
```
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
