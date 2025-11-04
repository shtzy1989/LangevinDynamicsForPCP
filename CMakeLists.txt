cmake_minimum_required(VERSION 3.16)
project(LangevinDynamics CXX)

include(CheckCXXCompilerFlag)

# ============================================================
# --- Require Intel MPI Compiler (mpiicpx) ---
# ============================================================
find_program(MPIICPX_EXECUTABLE mpiicpx)

if(NOT MPIICPX_EXECUTABLE)
    message(FATAL_ERROR
        "Intel MPI compiler 'mpiicpx' not found.\n"
        "Please load Intel MPI (e.g. source /opt/intel/oneapi/setvars.sh)\n"
        "or add mpiicpx to your PATH.")
endif()

set(CMAKE_CXX_COMPILER ${MPIICPX_EXECUTABLE})
message(STATUS "Using MPI compiler: ${CMAKE_CXX_COMPILER}")

# ============================================================
# --- Compiler settings ---
# ============================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ============================================================
# --- Detect AVX512 or AVX2 support ---
# ============================================================
check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512F)
check_cxx_compiler_flag("-mavx512dq" COMPILER_SUPPORTS_AVX512DQ)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)

if(COMPILER_SUPPORTS_AVX512F AND COMPILER_SUPPORTS_AVX512DQ)
    message(STATUS "AVX512 supported — using AVX512 optimization flags.")
    add_compile_options(-O2 -qopenmp -mavx2 -mavx512f -mavx512dq -fp-model=precise)
    set(SIMD_MODE "AVX512")
elseif(COMPILER_SUPPORTS_AVX2)
    message(STATUS "AVX512 not supported — using AVX2 optimization flags.")
    add_compile_options(-O2 -qopenmp -mavx2 -fp-model=precise -DSIMDLENGTH8)
    set(SIMD_MODE "AVX2")
else()
    message(FATAL_ERROR "Neither AVX512 nor AVX2 is supported by this compiler.")
endif()

# ============================================================
# --- Detect OpenMP support ---
# ============================================================
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM" OR
   CMAKE_CXX_COMPILER MATCHES "icpx" OR
   CMAKE_CXX_COMPILER MATCHES "mpiicpx")
    message(STATUS "Intel compiler detected → using -qopenmp")
    set(OPENMP_FLAGS "-qopenmp")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
       CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(STATUS "GCC/Clang detected → using -fopenmp")
    set(OPENMP_FLAGS "-fopenmp")
else()
    message(WARNING "Unknown compiler → OpenMP not enabled")
    set(OPENMP_FLAGS "")
endif()

set(CMAKE_C_FLAGS           "${CMAKE_C_FLAGS} ${OPENMP_FLAGS}")
set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} ${OPENMP_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${OPENMP_FLAGS}")

# ============================================================
# --- Include directories ---
# ============================================================
include_directories(${CMAKE_SOURCE_DIR}/src)

# ============================================================
# --- MKL and OpenMP Linking ---
# ============================================================
set(MKLROOT $ENV{MKLROOT})
if(NOT MKLROOT)
    message(FATAL_ERROR "MKLROOT environment variable not set. Please run: source /opt/intel/oneapi/setvars.sh")
endif()

set(MKL_LIB_DIR ${MKLROOT}/lib/intel64)
set(MKL_INC_DIR ${MKLROOT}/include)

include_directories(${MKL_INC_DIR})
link_directories(${MKL_LIB_DIR})

set(MKL_LIBRARIES
    mkl_blas95_ilp64
    mkl_intel_lp64
    mkl_core
    mkl_intel_thread
    mkl_gf_lp64
    mkl_rt
    mkl_gf_lp64
    mkl_gf_ilp64
)

# ============================================================
# --- Source files ---
# ============================================================
set(SRC_FILES
    src/vector.cpp
    src/DataLocator.cpp
    src/clsTime.cpp
    src/streamer.cpp
    src/argumentParser.cpp
    src/mpiKernel.cpp
    src/inputParser.cpp
    src/STATISTICS.cpp
    src/MOLUTILITY.cpp
    src/clsHistogram.cpp
    src/RANDOMNUMBER.cpp
    src/text.cpp
    src/simd.cpp
)

# ============================================================
# --- Libraries ---
# ============================================================
set(LINK_LIBS
    stdc++fs
    ${MKL_LIBRARIES}
)

# ============================================================
# --- Executables ---
# ============================================================
add_executable(langevinDynamicsSingleBarrier src/langevinDynamicsSingleBarrier.cpp ${SRC_FILES})
target_link_libraries(langevinDynamicsSingleBarrier PRIVATE ${LINK_LIBS})

add_executable(langevinDynamicsDoubleBarrier src/langevinDynamicsDoubleBarrier.cpp ${SRC_FILES})
target_link_libraries(langevinDynamicsDoubleBarrier PRIVATE ${LINK_LIBS})

# ============================================================
# --- Install rules ---
# ============================================================
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "Installation directory" FORCE)

install(TARGETS
    langevinDynamicsSingleBarrier
    langevinDynamicsDoubleBarrier
    RUNTIME DESTINATION bin
)

install(DIRECTORY src/
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")

# ============================================================
# --- Clean target ---
# ============================================================
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E rm -f *.o
            ${CMAKE_BINARY_DIR}/langevinDynamicsSingleBarrier
            ${CMAKE_BINARY_DIR}/langevinDynamicsDoubleBarrier
    COMMENT "Cleaning build files"
)

# ============================================================
# --- Summary ---
# ============================================================
message(STATUS "")
message(STATUS "============================================")
message(STATUS "  Build Configuration Summary")
message(STATUS "--------------------------------------------")
message(STATUS "  Compiler:   ${CMAKE_CXX_COMPILER}")
message(STATUS "  SIMD mode:  ${SIMD_MODE}")
message(STATUS "  MKL root:   ${MKLROOT}")
message(STATUS "  Libraries:  ${MKL_LIBRARIES}")
message(STATUS "  Install to: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "============================================")
message(STATUS "")

