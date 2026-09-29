# Tests, run with `ctest` (from the build directory, after building).
#
#   lattice_tessellation  Voronoi tessellation of a cubic lattice must be exact cubes.
#   synthetic_pipeline    end-to-end run on a synthetic SPH cloud; checks mass conservation and
#                         velocity recovery for sph-to-grid and the full Voronoi pipeline.
#                         Needs python3 (standard library only).

add_test(NAME lattice_tessellation COMMAND test_lattice_grid)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
  add_test(NAME synthetic_pipeline
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/synthetic_test.py
                   $<TARGET_FILE:kiwi> ${CMAKE_CURRENT_BINARY_DIR}/synthetic_work)
else()
  message(STATUS "python3 not found: skipping the synthetic_pipeline test")
endif()
