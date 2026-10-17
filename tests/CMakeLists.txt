# Tests, run with `ctest` (from the build directory, after building).
#
#   lattice_tessellation  Voronoi tessellation of a cubic lattice must be exact cubes.
#   synthetic_pipeline    end-to-end run on a synthetic SPH cloud; checks mass conservation and
#                         velocity recovery for sph-to-grid and the full Voronoi pipeline.
#                         Needs python3 (standard library only).
#   generic_kernel        generic-kernel analytic integral vs full_integral_new() (cubic) and brute-force quadrature.
#   kernel_spec_data_current  src/kernel_spec_data.cpp is what tools/kernel_tables/gen_cpp.py generates.
#   kernel_tables         derived kernel-integral coefficient tables vs brute-force quadrature
#                         (tools/kernel_tables). Needs python3 (standard library only).

add_test(NAME lattice_tessellation COMMAND test_lattice_grid)
add_test(NAME generic_kernel COMMAND test_generic_kernel)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
  add_test(NAME synthetic_pipeline
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/synthetic_test.py
                   $<TARGET_FILE:kiwi> ${CMAKE_CURRENT_BINARY_DIR}/synthetic_work)
  add_test(NAME kernel_tables
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/kernel_tables/verify.py 4)
  add_test(NAME kernel_spec_data_current
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/kernel_tables/gen_cpp.py --check)
else()
  message(STATUS "python3 not found: skipping the synthetic_pipeline and kernel_tables tests")
endif()
