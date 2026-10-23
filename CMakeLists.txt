cmake_minimum_required(VERSION 3.13)
project(kiwi CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

# --- OpenMP (mapping stage's particle loop) ---
# Optional, not required: the mapping stage runs single-threaded (correct,
# just slow) if OpenMP isn't found, since its #pragma omp lines are no-ops
# without it.
find_package(OpenMP)

# --- vendored voro++ ---
# Explicit list, matching voro++'s own src/Makefile `objs` exactly -- do NOT
# glob this directory: it also contains cmd_line.cc (voro++'s own CLI, has
# its own main()), v_base_wl.cc (a generated data table #included by
# v_base.cc, not a standalone TU), and voro++.cc (an all-in-one amalgamation
# that duplicates every symbol above if compiled alongside them).
set(VORO_SOURCES
  ${CMAKE_SOURCE_DIR}/third_party/voro/src/cell.cc
  ${CMAKE_SOURCE_DIR}/third_party/voro/src/common.cc
  ${CMAKE_SOURCE_DIR}/third_party/voro/src/container.cc
  ${CMAKE_SOURCE_DIR}/third_party/voro/src/unitcell.cc
  ${CMAKE_SOURCE_DIR}/third_party/voro/src/v_compute.cc
  ${CMAKE_SOURCE_DIR}/third_party/voro/src/c_loops.cc
  ${CMAKE_SOURCE_DIR}/third_party/voro/src/v_base.cc
  ${CMAKE_SOURCE_DIR}/third_party/voro/src/wall.cc
  ${CMAKE_SOURCE_DIR}/third_party/voro/src/pre_container.cc
  ${CMAKE_SOURCE_DIR}/third_party/voro/src/container_prd.cc
)
add_library(voro_static STATIC ${VORO_SOURCES})
target_include_directories(voro_static PUBLIC ${CMAKE_SOURCE_DIR}/third_party/voro/src)
# Vendored, unmodified: macOS flags the sprintf calls in voro++'s draw_pov() (POV-Ray export, never called
# by KIWI) as deprecated. Silence it for this target only rather than editing the upstream source.
target_compile_options(voro_static PRIVATE -Wno-deprecated-declarations)

# --- our library ---
file(GLOB KIWI_SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)
add_library(kiwi_core STATIC ${KIWI_SOURCES})
target_include_directories(kiwi_core PUBLIC
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/third_party
  ${CMAKE_SOURCE_DIR}/third_party/voro/src
)
target_link_libraries(kiwi_core PUBLIC voro_static)
# Lets diagnostics.cpp locate the bundled tools/plot_diagnostic.py by an
# absolute path regardless of the executable's cwd (e.g. running ./kiwi from
# inside build/, as the README's examples do).
target_compile_definitions(kiwi_core PUBLIC KIWI_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
if(OpenMP_CXX_FOUND)
  target_link_libraries(kiwi_core PUBLIC OpenMP::OpenMP_CXX)
else()
  message(WARNING "OpenMP not found -- the mapping stage will run single-threaded.")
endif()

# --- CLI executable ---
add_executable(kiwi src/main.cpp)
target_link_libraries(kiwi PRIVATE kiwi_core)
target_include_directories(kiwi PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/third_party)

# --- tests ---
enable_testing()
add_subdirectory(tests)

# --- permanent standalone tools (shipped, not part of the kiwi CLI itself) ---
add_executable(bench_grid_blocks tools/bench_grid_blocks.cpp)
target_link_libraries(bench_grid_blocks PRIVATE kiwi_core)
target_include_directories(bench_grid_blocks PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/third_party)

# Permanent regression test: Voronoi tessellation of a simple cubic lattice
# is exactly a partition into cubes -- see tools/test_lattice_grid.cpp. Also
# the ground-truth reference fixture for validating sph-to-grid (the direct
# SPH -> Cartesian mapping).
add_executable(test_lattice_grid tools/test_lattice_grid.cpp)
target_link_libraries(test_lattice_grid PRIVATE kiwi_core)
target_include_directories(test_lattice_grid PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/third_party)

# Regression test for the generic-kernel analytic integral (kiwi/kernel_integral_generic.hpp): the
# cubic spline must match the hand-written full_integral_new(), and all six kernels must match a
# brute-force quadrature of the kernel definition. See tools/test_generic_kernel.cpp.
add_executable(test_generic_kernel tools/test_generic_kernel.cpp)
target_link_libraries(test_generic_kernel PRIVATE kiwi_core)
target_include_directories(test_generic_kernel PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/third_party)

# Particle-only smoothing-length routine (Phantom-style neighbour-density
# iteration) -- first piece of the standalone inverse-mapping tool's particle
# machinery, built and validated independently before folding it in. See
# tools/h_from_neighbours.cpp's header comment. Links kiwi_core only for
# load_sph_snapshot() and OpenMP -- it does not touch voro++/tessellation at
# all.
add_executable(h_from_neighbours tools/h_from_neighbours.cpp)
target_link_libraries(h_from_neighbours PRIVATE kiwi_core)
target_include_directories(h_from_neighbours PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/third_party)

# Every shipped example config in configs/ must load -- see tools/test_configs.cpp.
add_executable(test_configs tools/test_configs.cpp)
target_link_libraries(test_configs PRIVATE kiwi_core)
target_include_directories(test_configs PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/third_party)
