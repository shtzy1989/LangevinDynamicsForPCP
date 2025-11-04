MPICXX   := mpiicpx

COMPILEFLAGS := -std=c++17 -O2 -mavx -mavx2 -fp-model=precise
COMPILEFLAGS_AVX512 := -mavx512f -mavx512dq -DSIMDLENGTH8

LIB          := -qopenmp -lstdc++fs -lmkl_intel_lp64 -lmkl_intel_lp64 \
          -lmkl_core \
          -lmkl_intel_thread \
          -lmkl_gf_lp64 \
          -lmkl_rt \
          -lmkl_gf_lp64 \
          -lmkl_gf_ilp64 

OBJ := vector.o \
       simd.o \
       text.o \
       DataLocator.o \
       clsTime.o \
       streamer.o \
       argumentParser.o \
       mpiKernel.o \
       inputParser.o \
       STATISTICS.o \
       MOLUTILITY.o \
       clsHistogram.o \
       RANDOMNUMBER.o

all:
	@echo "Please specify which target to build: 'make avx2' or 'make avx512'."
	@false

avx2: langevinDynamicsSingleBarrier_avx2 langevinDynamicsDoubleBarrier_avx2

%.o: %.cpp
	$(MPICXX) $(COMPILEFLAGS) -c $< -o $@

langevinDynamicsSingleBarrier_avx2: langevinDynamicsSingleBarrier.cpp $(OBJ)
	$(MPICXX) $(COMPILEFLAGS) -o $@ $^ $(LIB) 

langevinDynamicsDoubleBarrier_avx2: langevinDynamicsDoubleBarrier.cpp $(OBJ)
	$(MPICXX) $(COMPILEFLAGS) -o $@ $^ $(LIB) 

avx512: langevinDynamicsSingleBarrier_avx512 langevinDynamicsDoubleBarrier_avx512

%.o.avx512: %.cpp
	$(MPICXX) $(COMPILEFLAGS) $(COMPILEFLAGS_AVX512) -c $< -o $@

langevinDynamicsSingleBarrier_avx512: langevinDynamicsSingleBarrier.cpp $(OBJ:.o=.o.avx512)
	$(MPICXX) $(COMPILEFLAGS) $(COMPILEFLAGS_AVX512) -o $@ $^ $(LIB) 

langevinDynamicsDoubleBarrier_avx512: langevinDynamicsDoubleBarrier.cpp $(OBJ:.o=.o.avx512)
	$(MPICXX) $(COMPILEFLAGS) $(COMPILEFLAGS_AVX512) -o $@ $^ $(LIB) 

clean:
	rm -f *.o *.o.avx512 langevinDynamicsSingleBarrier langevinDynamicsDoubleBarrier

