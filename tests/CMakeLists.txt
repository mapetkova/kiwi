# Tests, run with `ctest` (from the build directory, after building).
#
#   lattice_tessellation  Voronoi tessellation of a cubic lattice must be exact cubes.
#   synthetic_pipeline    end-to-end run on a synthetic SPH cloud; checks mass conservation and
#                         velocity recovery for sph-to-grid and the full Voronoi pipeline.
#                         Needs python3 (standard library only).
#   generic_kernel        generic-kernel analytic integral vs full_integral_new() (cubic) and brute-force quadrature.
#   real_data             sph-to-grid on a real Phantom snapshot (tests/data/brick_03927.ascii, not in git, ~1.5 min):
#                         mass/momentum conservation for several kernels, generic == fast path. Skipped if absent.
#   real_data_cutout      the same on tests/data/brick_03927_cutout.ascii (x in [-36,-34], y in [-60,-58] pc, in git).
#   kernel_spec_data_current  src/kernel_spec_data.cpp is what tools/kernel_tables/gen_cpp.py generates.
#   phantom_dump          direct reader for Phantom binary dumps: synthetic dumps in many layouts vs the same
#                         particles as plain ascii; malformed dumps rejected (tests/phantom_dump_test.py).
#   example_configs       every configs/*.json loads through Config::load_from_file (tools/test_configs.cpp).
#   polaris_writer       `kiwi write-polaris` on a synthetic cloud, decoded with decode_polaris_grid.py: header,
#                         geometry (volumes, symmetric neighbours) and payload vs mapping.dat (polaris_writer_test.py).
#                         validate_writer_polaris.py, the one-off check against the legacy grid-file.cc output,
#                         needs that legacy file and is not run here.
#   bad_input             malformed configs and snapshots must give a one-line `[kiwi] error` and exit 1
#                         (bad_input_test.py).
#   kernel_tables         derived kernel-integral coefficient tables vs brute-force quadrature
#                         (tools/kernel_tables). Needs python3 (standard library only).

add_test(NAME lattice_tessellation COMMAND test_lattice_grid)
add_test(NAME generic_kernel COMMAND test_generic_kernel)

file(GLOB KIWI_EXAMPLE_CONFIGS ${CMAKE_SOURCE_DIR}/configs/*.json)
add_test(NAME example_configs COMMAND test_configs ${KIWI_EXAMPLE_CONFIGS})

find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
  add_test(NAME synthetic_pipeline
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/synthetic_test.py
                   $<TARGET_FILE:kiwi> ${CMAKE_CURRENT_BINARY_DIR}/synthetic_work)
  add_test(NAME phantom_dump
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/phantom_dump_test.py
                   $<TARGET_FILE:kiwi> ${CMAKE_CURRENT_BINARY_DIR}/phantom_dump_work)
  add_test(NAME polaris_writer
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/polaris_writer_test.py
                   $<TARGET_FILE:kiwi> ${CMAKE_CURRENT_BINARY_DIR}/polaris_writer_work)
  add_test(NAME bad_input
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/bad_input_test.py
                   $<TARGET_FILE:kiwi> ${CMAKE_CURRENT_BINARY_DIR}/bad_input_work)
  add_test(NAME kernel_tables
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/kernel_tables/verify.py 4)
  # Real simulation snapshot (tests/data, not tracked by git); skipped (exit 77) when it is not there.
  add_test(NAME real_data
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/real_data_test.py
                   $<TARGET_FILE:kiwi> ${CMAKE_CURRENT_SOURCE_DIR}/data/brick_03927.ascii
                   ${CMAKE_CURRENT_BINARY_DIR}/real_data_work)
  set_tests_properties(real_data PROPERTIES SKIP_RETURN_CODE 77 LABELS "slow" TIMEOUT 1800)
  # Same checks on a small spatial cutout of that snapshot (tracked in git, so it runs everywhere).
  add_test(NAME real_data_cutout
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/real_data_test.py
                   $<TARGET_FILE:kiwi> ${CMAKE_CURRENT_SOURCE_DIR}/data/brick_03927_cutout.ascii
                   ${CMAKE_CURRENT_BINARY_DIR}/real_data_cutout_work)
  set_tests_properties(real_data_cutout PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 600)
  add_test(NAME kernel_spec_data_current
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/kernel_tables/gen_cpp.py --check)
else()
  message(STATUS "python3 not found: skipping the synthetic_pipeline and kernel_tables tests")
endif()
