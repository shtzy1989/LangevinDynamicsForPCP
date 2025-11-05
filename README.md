# Langevin Dynamics for Porous Coordination Polymer (PCP)

This project provides two MPI-parallelized executables for Langevin dynamics simulations:
- **langevinDynamicsSingleBarrier** — langevin dynamics calculation of a single-barrier model  
- **langevinDynamicsDoubleBarrier** — langevin dynamics calculation of a double-barrier model  

---

## Directory Structure

```
src/              → all header and source files and Makefile
sample_input/     → sample input files
CMakeLists.txt    → main cmake build configuration file
install.sh        → quick script for installation
README.md         → README file
```

---

## Prerequisites

1. **Intel oneAPI**  
   Compiler: `mpiicpx`  
   Library: MAth Kernel Library (MKL) and Vector Statistical Library (VSL)  
   MPI: Intel MPI is recommended (`mpiicpx`)  
   Components required:  

   - Intel® oneAPI Base Toolkit  
   - Intel® oneAPI HPC Toolkit   
   
   Before building, set up the Intel oneAPI environment:
   ```bash
   source {path_to_oneapi}/setvars.sh
   ```

2. **CMake ≥ 3.16**

3. **CPUs with AVX512 or AVX2 support**  
For performance reasons, plain C++ version has not been implemented  
---

## Build Instructions

### Option 1: Build using CMake (recommended)

The provided `CMakeLists.txt` automatically:
- Detects the compiler (`mpiicpx`)
- Detects supported SIMD instructions (AVX2 / AVX512)
- Detects MKL and VSL

Steps:
```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=install
cd build
make -j
make install
```

After installation, the executables are found by default in:
```
install/bin/
```

---

### Option 2: Build using Makefiles

If you prefer, the makefile is available under `src/`:
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
Thus generated executable files will be named with suffix, `_avx512` or `_avx2`.

---

## Run Instructions

Run the executables to generate default input files:
```bash
./langevinDynamicsSingleBarrier -o input_single.inp
./langevinDynamicsSingleBarrier -o input_double.inp
```
This command can be used to correct options in existing input files or reformat them.  
The existing input file will be backed up by renaming.  
Refer to sample input files under `sample_input` for examples.  

Then, use `mpirun` or `mpiexec` to launch either executable:
```bash
mpirun -np 4 ./langevinDynamicsSingleBarrier -i input_single.inp
mpirun -np 4 ./langevinDynamicsDoubleBarrier -i input_double.inp
```
Prediction of execution time will be displayed if PRINTTIMER is set to true.

---

## Physical Models  

The underdamped langevin equation reads,  
$m \frac{d^2\mathbf{x}}{dt^2} = -(\frac{dV(\mathbf{x})}{dx}+\frac{d \Delta V(\mathbf{x})}{dx}\eta(t)) - \gamma_{s} m \frac{dx}{dt} + \sqrt{2 \gamma_s m k_B T} \, \boldsymbol{\xi}(t).$

The overdamped langevin equation reads,  
$m \gamma_s \frac{d\mathbf{x}}{dt} = -(\frac{dV(\mathbf{x})}{dx}+\frac{d \Delta V(\mathbf{x})}{dx}\eta(t)) + \sqrt{2 \gamma_s m k_B T} \, \boldsymbol{\xi}(t).$  

The dichotomic noise $\eta(t)$ is governed by,  
$\langle \eta(t) \rangle=0$, and  
$\langle \eta(t) \eta(s) \rangle=\mathrm{exp}(-2\gamma_\eta|t-s|)$.  

---

## Input Options  
- **OUTPUTBASE**  
The base of the output file names.
- **OVERDAMPED**  
true = use overdamped langevin equation.  
false = use underdamped langevin equation.
- **X0**  
Initial position of particles.
- **V0**  
Initial velocity of particles (irrelevant when using overdamped langevin equation).  
- **GAMMA**  
The friction coefficient $\gamma_s$.  
- **KB**  
The Boltzmann coefficient $k_B$.  
- **TEMPERATURE**  
The tempeature $T$.  
- **MASS**  
The mass of the particle $m$.  
- **RANDOMIZEV**  
true = randomize the initial velocity of particles, overrides V0.  
false = not to randomize the initial velocity of particles and use V0.  
- **PERIODICLENGTH**  
The periodic length $L$ of the potential energy.  
The primitive period is located between $[-L/2, L/2)$  
- **BARRIERHEIGHTS**  
Heights of the barriers in the two states, between which the fluctuation occurs.  
Requires two values.  
Only legal in langevinDynamicsSingleBarrier.  
- **BARRIERHEIGHTS0**  
Heights of the first barriers in the two states, between which the fluctuation occurs.  
Requires two values, with the first higher than the second.  
Only legal in langevinDynamicsDoubleBarrier.  
- **BARRIERHEIGHTS1**  
Heights of the second barriers in the two states, between which the fluctuation occurs.  
Requires two values, with the first higher than the second.  
Only legal in langevinDynamicsDoubleBarrier.  
- **BARRIERHALFWIDTH**  
Half width of the barrier.  
Only legal in langevinDynamicsSingleBarrier.  
- **BARRIERHALFWIDTH0**  
Half width of the first barrier.  
Only legal in langevinDynamicsDoubleBarrier.  
- **BARRIERHALFWIDTH1**  
Half width of the second barrier.  
Only legal in langevinDynamicsDoubleBarrier.  
- **BARRIERPOSITION0**  
Position of the first barrier in the primitiva period.  
Only legal in langevinDynamicsDoubleBarrier.  
The barrier in the single barrier model is always positioned at $x = nL, n \in \mathbb{Z}$
- **BARRIERPOSITION1**  
Position of the second barrier in the primitiva period.  
Only legal in langevinDynamicsDoubleBarrier.  
- **TIMESTEP**  
Timestep used in the integration.  
- **NUMBEROFSTEPS**  
The number of steps of each run.  
- **NUMBEROFPREVIOUSSTEPS**  
The number of previous steps for performing restarting runs.  
- **NUMBEROFRUNS**  
The number of consecutive runs bundled in a single calculation.  
Note that the number of parallel runs is determined by number of MPI threads.  
- **FREQSAVE**  
The frequency for saving trajectories.  
It determines the minimal timescale of the calculated MSD.  
- **FREQSTAT**  
The frquency for first passage time, transition path time, and waiting time calculations.  
- **DISTRIBUTIONV0**  
true = calculate distribution of initial velocity (irrelevant when using overdamped langevin equation).  
false = not to calculate distribution of initial velocity.  
- **DISTRIBUTIONS0**  
true = calculate distribution of initial barrier state.  
false = not to calculate distribution of initial barrier state.  
- **TAU**  
Timescale $\tau_\eta=1/\gamma_\eta$ of the bichotomic noise $\eta(t)$.
- **INITIALSTATE**  
Initial barrier state (0 or 1, irrelevant when randomizing the initial state).  
- **RANDOMIZESTATES**  
true = randomize the initial barrier state.  
false = not to randomize the initial barrier state.  
- **UNIFORMINITALSTATES**  
true = use the initial barrier state at x0 for all periodic barriers.  
false = randomize the periodic barriers other than the one at x0.  
- **EQ**  
The equilibrium constant from barrier state 0 to state 1.  
- **WAITINGTIME**  
true = calculate the waiting time between barrier crossing events.  
false = not to calculate the waiting time between barrier crossing events.  
- **FPT**  
true = calculate the first passage time.  
false = not to calculate the first passage time.  
- **TPT**  
true = calculate the transition path time.  
false = not to calculate the transition path time.  
- **FPTCG**  
true $N$ = calculate the first passage time between $N$ barriers.  
false $N$ = not to do the calculation.  
Can be specified multiple times for different $N$ values.  
- **READRESTART**  
true = read a restart file.  
false = not to read a restart file.  
- **WRITERESTART**  
true = write a restart file at the end of the calculation.  
false = not to write a restart file at the end of the calculation.  
- **INPUTRESTART**  
Input restart file name base (corresponding to the part in OUTPUTBASE).  
- **OUTPUTRESTART**  
Output restart file name base (corresponding to the part in OUTPUTBASE).  
- **WRITEX**  
true = output trajectories.  
false = not to output trajectories.  
- **NUMBEROFLOGSAMPLES**  
The number of sample points per order in the log scale MSD log file.  
- **RANDOMSEED**  
The master random number seed used to generate the random number seed of each mpi thread.  
- **PRINTTIMER**  
true = print the time prediction during calculation.  
false = not to print the time prediction during calculation.  

---

## Contact

For build or usage questions, please email the author with the follow information:
- Compiler and MKL version  
- CMake output from configuration (`cmake ..`)  
- Any build log messages

---