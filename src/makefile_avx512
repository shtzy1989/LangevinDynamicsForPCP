MPICXX   := mpiicpx

# AVX2
#COMPILEFLAGS := -std=c++17 -O2 -mavx2 -fp-model=precise -DSIMDLENGTH8 

# AVX512
COMPILEFLAGS := -std=c++17 -O2 -mavx2 -mavx512f -mavx512dq -fp-model=precise #-DSIMDLENGTH8 

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

all: langevinDynamicsSingleBarrier langevinDynamicsDoubleBarrier

%.o: %.cpp
	$(MPICXX) $(COMPILEFLAGS) -c $< -o $@

langevinDynamicsSingleBarrier: langevinDynamicsSingleBarrier.cpp $(OBJ)
	$(MPICXX) $(COMPILEFLAGS) -o $@ $^ $(LIB) 

langevinDynamicsDoubleBarrier: langevinDynamicsDoubleBarrier.cpp $(OBJ)
	$(MPICXX) $(COMPILEFLAGS) -o $@ $^ $(LIB) 

clean:
	rm -f *.o langevinDynamicsSingleBarrier langevinDynamicsDoubleBarrier

