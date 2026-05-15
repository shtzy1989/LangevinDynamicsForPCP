# Langevin Dynamics for Porous Coordination Polymer (PCP)

This project provides two MPI-parallelized executables for Langevin dynamics simulations:
- **langevinDynamicsSingleBarrier** — Langevin dynamics calculation of a single-barrier model  
- **langevinDynamicsDoubleBarrier** — Langevin dynamics calculation of a double-barrier model  

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
   Library: Math Kernel Library (MKL) and Vector Statistical Library (VSL)  
   MPI: Intel MPI is recommended (`mpiicpx` and `mpirun`)  
   Components required:  

   - Intel® oneAPI Base Toolkit  
   - Intel® oneAPI HPC Toolkit   
   
   Before building, set up the Intel oneAPI environment:
   ```bash
   source {path_to_oneapi}/setvars.sh
   ```

2. **CMake ≥ 3.16**

3. **CPUs with AVX512 or AVX2 support**  
For performance reasons, a plain C++ version has not been implemented  
---

## Build Instructions

Installation typically takes only a few minutes.

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

## Sample Inputs

Example input files are provided under the `sample_input` directory.

Folder structure:
```text
sample_input/
├── quick_test/
│   ├── 2.0_7.2_m1.0_1.00e-05.inp
│   ├── run.sh
│   └── executed_job/
│       └── ...(see Expected Outputs)
├── CuIDB_double/
│   ├── 2.0_7.2_m1.0_1.00e-05.inp
│   ├── 3.0_10.8_m1.0_1.00e-05.inp
│   ├── 4.0_14.4_m1.0_1.00e-05.inp
│   └── 5.0_18.0_m1.0_1.00e-05.inp
├── CuIDB_single/
│   ├── 2.0_7.2_m1.0_1.00e-05.inp
│   ├── 3.0_10.8_m1.0_1.00e-05.inp
│   ├── 4.0_14.4_m1.0_1.00e-05.inp
│   └── 5.0_18.0_m1.0_1.00e-05.inp
└── CuOPTz/
    ├── 3.50_9.45_m1.0_1.00e-05.inp
    ├── 4.50_12.15_m1.0_1.00e-05.inp
    ├── 5.50_14.85_m1.0_1.00e-05.inp
    └── 6.50_17.55_m1.0_1.00e-05.inp
```

The `quick_test` directory contains a quick test calculation intended to verify that the program is installed and functioning correctly. The `executed_job` subdirectory contains the output files produced by performing the calculation using the run.sh.

The `CuIDB_single` and `CuOPTz` directories contain input files for calculations using `langevinDynamicsSingleBarrier`.

The `CuIDB_double` directory contains input files for calculations using `langevinDynamicsDoubleBarrier`.

---

## Run Instructions

### Generating Input Files

Run the executables to generate default input files:
```bash
./langevinDynamicsSingleBarrier -o input_single.inp
./langevinDynamicsDoubleBarrier -o input_double.inp
```

These commands can also be used to:
- update missing or obsolete input options,
- reformat existing input files.

Existing input files are automatically backed up by renaming.

Refer to the example input files under `sample_input` for reference.

---

### Running Calculations

Use `mpirun` or `mpiexec` to launch the simulations:
```bash
mpirun -np 4 ./langevinDynamicsSingleBarrier -i input_single.inp

mpirun -np 4 ./langevinDynamicsDoubleBarrier -i input_double.inp
```

An estimated execution time is displayed if `PRINTTIMER` is set to `true` in the input file.

---

### Running the Quick Test

After installation, begin with the quick test under `sample_input/quick_test`.

Run the test calculation inside the `quick_test` directory:
```bash
bash ./run.sh
```

The script launches four MPI calculations, each using either four or eight OpenMP threads depending on whether AVX2 or AVX512 SIMD instructions are employed, resulting in a total of 16 or 32 independent simulations.

---

### Reference Test Environments

#### Test Environment 1

- CPU: Intel(R) Xeon(R) Gold 6136 CPU @ 3.00GHz
- OS: Ubuntu 24.04 LTS
- Compiler: Intel(R) oneAPI DPC++/C++ Compiler 2024.1.0 (2024.1.0.20240308)
- MPI: Intel(R) MPI Library for Linux* OS, Version 2021.12 Build 20240213
- SIMD: AVX + SVML
- Wall-clock time: 7.2 min

#### Test Environment 2

- CPU: Intel(R) Core(TM) Ultra 7 265U
- OS: Ubuntu 24.04 LTS (WSL on Windows 11)
- Compiler: Intel(R) oneAPI DPC++/C++ Compiler 2026.0.0 (2026.0.0.20260331)
- MPI: Intel(R) MPI Library for Linux* OS, Version 2021.18.0 Build 20260327
- SIMD: AVX + SVML
- Wall-clock time: 5.7 min

---

### Production Calculations

For production calculations:
- Use `langevinDynamicsSingleBarrier` for input files under `CuIDB_single` and `CuOPTz`.
- Use `langevinDynamicsDoubleBarrier` for input files under `CuIDB_double`.

Example:
```bash
mpirun -np 4 ./langevinDynamicsSingleBarrier -i ./2.0_7.2_m1.0_1.00e-05.inp

mpirun -np 4 ./langevinDynamicsDoubleBarrier -i ./2.0_7.2_m1.0_1.00e-05.inp
```
---

## Physical Models  

The underdamped Langevin equation reads,  
$m \frac{d^2\mathbf{x}}{dt^2} = -(\frac{dV(\mathbf{x})}{dx}+\frac{d \Delta V(\mathbf{x})}{dx}\eta(t)) - \gamma_{s} m \frac{dx}{dt} + \sqrt{2 \gamma_s m k_B T} \, \boldsymbol{\xi}(t).$

The overdamped Langevin equation reads,  
$m \gamma_s \frac{d\mathbf{x}}{dt} = -(\frac{dV(\mathbf{x})}{dx}+\frac{d \Delta V(\mathbf{x})}{dx}\eta(t)) + \sqrt{2 \gamma_s m k_B T} \, \boldsymbol{\xi}(t).$  

The dichotomic noise $\eta(t)$ is governed by,  
$\langle \eta(t) \rangle=0$, and  
$\langle \eta(t) \eta(s) \rangle=\mathrm{exp}(-2\gamma_\eta|t-s|)$.  

---

## Langevin Dynamics Algorithms

Both underdamped and overdamped Langevin dynamics simulations were performed in a one-dimensional periodically repeating two-state barrier potential. The barrier state switched stochastically between two barrier heights according to prescribed transition probabilities. The overdamped version is used primarily.

### Underdamped Langevin Dynamics

In the underdamped case, both particle positions and velocities were propagated using a velocity-Verlet-type Langevin integrator:

```text
Initialize x, v
Initialize barrier state
Compute initial force F(x, barrier)

For each timestep:

    # Half-step velocity update
    v ← v + (F / m) * dt/2

    # Half-step position update
    x ← x + v * dt/2

    # Langevin thermostat
    R ← Gaussian random number
    v ← exp(-γ dt) * v
         + sqrt(kB T (1 - exp(-2γdt)) / m) * R

    # Second half-step position update
    x ← x + v * dt/2

    # Update periodic well index
    Determine current periodic cell and well

    # Barrier-state transition
    Switch barrier state according to transition probabilities

    # Update barrier parameters
    Update barrier height and slopes

    # Compute new deterministic force
    Compute F(x, barrier)

    # Final half-step velocity update
    v ← v + (F / m) * dt/2

    # Record observables
    Update statistics and save trajectory if required
```

### Overdamped Langevin Dynamics

In the overdamped limit, inertial effects were neglected and only particle positions were propagated:

```text
Initialize x
Initialize barrier state
Compute initial force F(x, barrier)

For each timestep:

    # Thermal fluctuation
    R ← Gaussian random number

    # Overdamped Langevin update
    x ← x
         + (F / (mγ)) * dt
         + sqrt(2kBTdt / (mγ)) * R

    # Update periodic well index
    Determine current periodic cell and well

    # Barrier-state transition
    Switch barrier state according to transition probabilities

    # Update barrier parameters
    Update barrier height and slopes

    # Compute new deterministic force
    Compute F(x, barrier)

    # Record observables
    Update statistics and save trajectory if required
```

The periodic barrier potential consisted of a piecewise linear two-state profile. At each timestep, the barrier state could switch between two barrier heights according to stochastic transition probabilities determined by the characteristic switching timescale. During the simulations, quantities such as trajectories, waiting-time distributions, first-passage times, and transition-path statistics were accumulated for analysis.

---

## Input Options  
- **OUTPUTBASE**  
The base of the output file names.
- **OVERDAMPED**  
true = use overdamped Langevin equation.  
false = use underdamped Langevin equation.
- **X0**  
Initial position of particles.
- **V0**  
Initial velocity of particles (irrelevant when using overdamped Langevin equation).  
- **GAMMA**  
The friction coefficient $\gamma_s$.  
- **KB**  
The Boltzmann coefficient $k_B$.  
- **TEMPERATURE**  
The temperature $T$.  
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
Position of the first barrier in the primitive period.  
Only legal in langevinDynamicsDoubleBarrier.  
The barrier in the single barrier model is always positioned at $x = nL, n \in \mathbb{Z}$
- **BARRIERPOSITION1**  
Position of the second barrier in the primitive period.  
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
The frequency for first passage time, transition path time, and waiting time calculations.  
- **DISTRIBUTIONV0**  
true = calculate distribution of initial velocity (irrelevant when using overdamped Langevin equation).  
false = not to calculate distribution of initial velocity.  
- **DISTRIBUTIONS0**  
true = calculate distribution of initial barrier state.  
false = not to calculate distribution of initial barrier state.  
- **TAU**  
Timescale $\tau_\eta=1/\gamma_\eta$ of the dichotomic noise $\eta(t)$.
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
true $N$ = calculate the coarse-grained first passage time across $N$ barriers.  
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

## Expected Outputs

### Dynamical Quantities and Timing Information

```text
m1.0_t1.00e-05_0_msd.log
m1.0_t1.00e-05_0_waitingTime.log
m1.0_t1.00e-05_0_time.log
```

- `msd.log`: Main log file for the calculated mean squared displacement.
- `waitingTime.log`: Histogram of waiting times between barrier crossing events.
- `time.log`: Execution time and performance information.


### Restart Files

```text
m1.0_t1.00e-05_restart_0_00000_00000.restart
m1.0_t1.00e-05_restart_0_00001_00000.restart
m1.0_t1.00e-05_restart_0_00002_00000.restart
m1.0_t1.00e-05_restart_0_00003_00000.restart
```

- Binary restart files generated for individual MPI ranks.
- These files can be used to resume interrupted simulations.
- When restarting a calculation, put the part corresponding to m1.0_t1.00e-05_restart_0 in INPUTRESTART option in the input.


### First Passage Time: Periodic Boundary to Nearest Midpoint

```text
m1.0_t1.00e-05_0_firstPassageTime.dat
m1.0_t1.00e-05_0_firstPassageTime.log
```

Binary data and log file for the first passage time from a periodic boundary to the nearest midpoint within the period.

### First Passage Time: Periodic Boundary to Neighboring Boundary

```text
m1.0_t1.00e-05_0_firstPassageTime2.dat
m1.0_t1.00e-05_0_firstPassageTime2.log
```

Binary data and log file for the first passage time from a periodic boundary to the neighboring periodic boundary.

### Coarse-Grained First Passage Time

```text
m1.0_t1.00e-05_0_firstPassageTimeCG2_0000_0002.dat
m1.0_t1.00e-05_0_firstPassageTimeCG2_0001_0003.dat
m1.0_t1.00e-05_0_firstPassageTimeCG2_0002_0004.dat
m1.0_t1.00e-05_0_firstPassageTimeCG2_0003_0005.dat
m1.0_t1.00e-05_0_firstPassageTimeCG2_0004_0010.dat
m1.0_t1.00e-05_0_firstPassageTimeCG2_0005_0100.dat
m1.0_t1.00e-05_0_firstPassageTimeCG2_0006_1000.dat
m1.0_t1.00e-05_0_firstPassageTimeCG2.log
```

- `*.dat`: Binary data for coarse-grained first passage times from a periodic boundary to distant neighboring boundary.
- The suffixes indicate the target neighbor distance (e.g., 2nd, 3rd, 10th, 100th, and 1000th neighbors).
- `*.log`: Log file for the coarse-grained first passage time analysis.

### Transition Path Time: Periodic Boundary to Nearest Midpoint

```text
m1.0_t1.00e-05_0_transitionPathTime.dat
m1.0_t1.00e-05_0_transitionPathTime.log
```

Binary data and log file for the transition path time from a periodic boundary to the nearest midpoint within the period.

### Transition Path Time: Periodic Boundary to Neighboring Boundary

```text
m1.0_t1.00e-05_0_transitionPathTime2.dat
m1.0_t1.00e-05_0_transitionPathTime2.log
```

Binary data and log file for the transition path time from a periodic boundary to the neighboring periodic boundary.

### Barrier Height and State Distributions

```text
m1.0_t1.00e-05_0_initialBarrierHeightHistogram.log
m1.0_t1.00e-05_0_instantBarrierHeightHistogram.log
m1.0_t1.00e-05_0_State0Distribution.log
m1.0_t1.00e-05_0_stateCount.log
```

- `initialBarrierHeightHistogram.log`: Histogram of barrier heights when the particle is located at a periodic boundary.
- `instantBarrierHeightHistogram.log`: Histogram of barrier heights at the moment of barrier crossing.
- `State0Distribution.log`: Histogram of the initial barrier states.
- `stateCount.log`: Statistics for the number and fraction of each barrier state.

### Output Trajectories (optional, when WRITEX is true)
```text
m1.0_t1.00e-05_0_trajectory_0000000000_0000000000_0000000000.log
m1.0_t1.00e-05_0_trajectory_0000000000_0000000000_0000000001.log
...
```

---

## Contact

For build or usage questions, please email the author with the following information:
- Compiler and MKL version  
- CMake output from configuration (`cmake ..`)  
- Any build log messages

---
