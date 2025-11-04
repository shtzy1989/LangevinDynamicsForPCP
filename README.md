# Langevin Dynamics Simulation

This project provides two MPI-parallelized executables for Langevin dynamics simulations:
- **langevinDynamicsSingleBarrier** — simulation of a single-barrier model  
- **langevinDynamicsDoubleBarrier** — simulation of a double-barrier model  

Both programs use the Vector Statistical Library (VSL) from Intel MKL for random number generation.

---

## Directory Structure

```
src/              → all source and header files
sample_input/     → sample input files of single/double barrier [Cu(IDB)] models, and [Cu(OPTz)] model
CMakeLists.txt    → main build configuration file
```

---

## Prerequisites

1. **Intel oneAPI** (providing `mpiicpx`, MKL, and VSL)  
   Before building, set up the Intel environment:
   ```bash
   source /opt/intel/oneapi/setvars.sh
   ```

2. **CMake ≥ 3.16**

3. **MPI environment**  
   Intel MPI is recommended; `mpiicpx` will be automatically detected.

4. CPUs with AVX512 or AVX2 support (for performance reasons, plain C++ version has not been implemented)
---

## Build Instructions

### Option 1: Build using CMake (recommended)

The provided `CMakeLists.txt` automatically:
- Detects your compiler (`mpiicpx`, `icpx`, GCC, or Clang`)
- Detects supported SIMD instructions (AVX2 / AVX512) automatically
- Links MKL and VSL correctly without manual flags

Steps:
```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=install
cd build
make -j
make install
```

After installation, the executables are found in:
```
install/bin/
```

---

### Option 2: Build using Makefiles (manual control)

If you prefer, the original Makefiles are available under `src/`:
- `makefile`

To build manually with avx2:
```bash
cd src
make avx2
```
To build manually with avx512:
```bash
cd src
make avx512
```
The Makefile uses `mpiicpx` and MKL flags explicitly and provide full control of compiler options.

---

## Run Instructions

Run the executables to generate default input files:
```bash
./langevinDynamicsSingleBarrier -o input_single.inp
./langevinDynamicsSingleBarrier -o input_double.inp
```
Use `mpirun` or `mpiexec` to launch either executable:
```bash
mpirun -np 4 ./langevinDynamicsSingleBarrier -i input_single.inp
mpirun -np 4 ./langevinDynamicsDoubleBarrier -i input_double.inp
```
Each program reads its own input file for system parameters.

---

## Notes

- The current CMakeLists.txt links all required MKL and VSL components  
  (`mkl_core`, `mkl_intel_lp64`, `mkl_intel_thread`, `mkl_vml_def`, `mkl_vsl_core`, `mkl_vsl_def`, and `iomp5`) automatically.
- The user does **not** need to edit compiler flags or specify AVX2/AVX512 options manually.
- Installation path defaults to `<project_root>/install`, and can be overridden using:
  ```bash
  cmake -DCMAKE_INSTALL_PREFIX=/your/path ..
  ```

---

## Contact

For build or usage questions, please include:
- Compiler and MKL version  
- CMake output from configuration (`cmake ..`)  
- Any build log messages

