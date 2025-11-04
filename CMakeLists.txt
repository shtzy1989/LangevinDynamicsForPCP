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

include(CheckCXXCompilerFlag)
include(CheckCXXSourceRuns)

set(CMAKE_REQUIRED_FLAGS "-O2")

# Step 1: check if compiler *accepts* the flags
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512F)

# Step 2: test at runtime if the CPU supports these instruction sets
if(COMPILER_SUPPORTS_AVX512F)
    message(STATUS "Compiler supports AVX512F, checking CPU capability...")

    file(WRITE ${CMAKE_BINARY_DIR}/check_avx512.cpp
        "#include <immintrin.h>\n"
        "int main() {\n"
        "#ifdef __GNUC__\n"
        "  if (!__builtin_cpu_supports(\"avx512f\")) return 1;\n"
        "#endif\n"
        "  return 0;\n"
        "}\n"
    )
    try_run(RUN_RESULT COMPILE_RESULT
        ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}/check_avx512.cpp
        COMPILE_DEFINITIONS "-mavx512f"
    )

    if(COMPILE_RESULT AND RUN_RESULT EQUAL 0)
        set(USE_AVX512 TRUE)
        message(STATUS "CPU supports AVX512F — enabling -mavx512f")
    else()
        message(STATUS "CPU does not support AVX512F — skipping")
    endif()
endif()

if(NOT USE_AVX512 AND COMPILER_SUPPORTS_AVX2)
    message(STATUS "Compiler supports AVX2, checking CPU capability...")

    file(WRITE ${CMAKE_BINARY_DIR}/check_avx2.cpp
        "#include <immintrin.h>\n"
        "int main() {\n"
        "#ifdef __GNUC__\n"
        "  if (!__builtin_cpu_supports(\"avx2\")) return 1;\n"
        "#endif\n"
        "  return 0;\n"
        "}\n"
    )
    try_run(RUN_RESULT2 COMPILE_RESULT2
        ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}/check_avx2.cpp
        COMPILE_DEFINITIONS "-mavx2"
    )

    if(COMPILE_RESULT2 AND RUN_RESULT2 EQUAL 0)
        set(USE_AVX2 TRUE)
        message(STATUS "CPU supports AVX2 — enabling -mavx2")
    else()
        message(STATUS "CPU does not support AVX2 — skipping")
    endif()
endif()

# Step 3: set appropriate compiler flags
if(USE_AVX512)
    set(AVX_FLAGS "-mavx512f -mavx512dq")
elseif(USE_AVX2)
    set(AVX_FLAGS "-mavx2")
else()
    message(FATAL_ERROR "❌ Neither AVX2 nor AVX512F are supported on this machine.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${AVX_FLAGS}")
message(STATUS "Final SIMD flags: ${AVX_FLAGS}")

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

set(CMAKE_C_FLAGS           "${CMAKE_C_FLAGS} ${OPENMP_FLAGS} -fp-model=precise")
set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} ${OPENMP_FLAGS} -fp-model=precise")
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

