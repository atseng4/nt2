# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.6

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jt/DevC++/dev_lasmea/github/nt2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jt/DevC++/dev_lasmea/github/nt2

# Include any dependencies generated for this target.
include include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/depend.make

# Include the progress variables for this target.
include include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/progress.make

# Include the compile flags for this target's objects.
include include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/flags.make

include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o: include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/flags.make
include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o: include/nt2/toolbox/gsl_specfun/bench/simd/gsl_sf_bessel_i0.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jt/DevC++/dev_lasmea/github/nt2/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o"
	cd /home/jt/DevC++/dev_lasmea/github/nt2/include/nt2/toolbox/gsl_specfun/bench/simd && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS)  -mssse3 -msse3 -msse2 -mmmx -o CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o -c /home/jt/DevC++/dev_lasmea/github/nt2/include/nt2/toolbox/gsl_specfun/bench/simd/gsl_sf_bessel_i0.cpp

include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.i"
	cd /home/jt/DevC++/dev_lasmea/github/nt2/include/nt2/toolbox/gsl_specfun/bench/simd && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS)  -mssse3 -msse3 -msse2 -mmmx -E /home/jt/DevC++/dev_lasmea/github/nt2/include/nt2/toolbox/gsl_specfun/bench/simd/gsl_sf_bessel_i0.cpp > CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.i

include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.s"
	cd /home/jt/DevC++/dev_lasmea/github/nt2/include/nt2/toolbox/gsl_specfun/bench/simd && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS)  -mssse3 -msse3 -msse2 -mmmx -S /home/jt/DevC++/dev_lasmea/github/nt2/include/nt2/toolbox/gsl_specfun/bench/simd/gsl_sf_bessel_i0.cpp -o CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.s

include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o.requires:
.PHONY : include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o.requires

include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o.provides: include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o.requires
	$(MAKE) -f include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/build.make include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o.provides.build
.PHONY : include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o.provides

include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o.provides.build: include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o
.PHONY : include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o.provides.build

# Object files for target gsl_sf_bessel_i0.gsl_specfun.simd.bench
gsl_sf_bessel_i0_gsl_specfun_simd_bench_OBJECTS = \
"CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o"

# External object files for target gsl_sf_bessel_i0.gsl_specfun.simd.bench
gsl_sf_bessel_i0_gsl_specfun_simd_bench_EXTERNAL_OBJECTS =

include/nt2/toolbox/gsl_specfun/bench/simd/gsl_sf_bessel_i0.gsl_specfun.simd.bench: include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o
include/nt2/toolbox/gsl_specfun/bench/simd/gsl_sf_bessel_i0.gsl_specfun.simd.bench: lib/Release/libnt2.so.3.0.0
include/nt2/toolbox/gsl_specfun/bench/simd/gsl_sf_bessel_i0.gsl_specfun.simd.bench: include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/build.make
include/nt2/toolbox/gsl_specfun/bench/simd/gsl_sf_bessel_i0.gsl_specfun.simd.bench: include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable gsl_sf_bessel_i0.gsl_specfun.simd.bench"
	cd /home/jt/DevC++/dev_lasmea/github/nt2/include/nt2/toolbox/gsl_specfun/bench/simd && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/build: include/nt2/toolbox/gsl_specfun/bench/simd/gsl_sf_bessel_i0.gsl_specfun.simd.bench
.PHONY : include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/build

include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/requires: include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/gsl_sf_bessel_i0.cpp.o.requires
.PHONY : include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/requires

include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/clean:
	cd /home/jt/DevC++/dev_lasmea/github/nt2/include/nt2/toolbox/gsl_specfun/bench/simd && $(CMAKE_COMMAND) -P CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/cmake_clean.cmake
.PHONY : include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/clean

include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/depend:
	cd /home/jt/DevC++/dev_lasmea/github/nt2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jt/DevC++/dev_lasmea/github/nt2 /home/jt/DevC++/dev_lasmea/github/nt2/include/nt2/toolbox/gsl_specfun/bench/simd /home/jt/DevC++/dev_lasmea/github/nt2 /home/jt/DevC++/dev_lasmea/github/nt2/include/nt2/toolbox/gsl_specfun/bench/simd /home/jt/DevC++/dev_lasmea/github/nt2/include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : include/nt2/toolbox/gsl_specfun/bench/simd/CMakeFiles/gsl_sf_bessel_i0.gsl_specfun.simd.bench.dir/depend

