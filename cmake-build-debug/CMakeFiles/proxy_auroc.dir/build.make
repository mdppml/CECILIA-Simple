# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/meteakgun/CLionProjects/CECILIA

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/meteakgun/CLionProjects/CECILIA/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/proxy_auroc.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/proxy_auroc.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/proxy_auroc.dir/flags.make

CMakeFiles/proxy_auroc.dir/apps/auroc/proxy.cpp.o: CMakeFiles/proxy_auroc.dir/flags.make
CMakeFiles/proxy_auroc.dir/apps/auroc/proxy.cpp.o: ../apps/auroc/proxy.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/meteakgun/CLionProjects/CECILIA/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/proxy_auroc.dir/apps/auroc/proxy.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/proxy_auroc.dir/apps/auroc/proxy.cpp.o -c /Users/meteakgun/CLionProjects/CECILIA/apps/auroc/proxy.cpp

CMakeFiles/proxy_auroc.dir/apps/auroc/proxy.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/proxy_auroc.dir/apps/auroc/proxy.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/meteakgun/CLionProjects/CECILIA/apps/auroc/proxy.cpp > CMakeFiles/proxy_auroc.dir/apps/auroc/proxy.cpp.i

CMakeFiles/proxy_auroc.dir/apps/auroc/proxy.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/proxy_auroc.dir/apps/auroc/proxy.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/meteakgun/CLionProjects/CECILIA/apps/auroc/proxy.cpp -o CMakeFiles/proxy_auroc.dir/apps/auroc/proxy.cpp.s

CMakeFiles/proxy_auroc.dir/utils/parse_options.cpp.o: CMakeFiles/proxy_auroc.dir/flags.make
CMakeFiles/proxy_auroc.dir/utils/parse_options.cpp.o: ../utils/parse_options.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/meteakgun/CLionProjects/CECILIA/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/proxy_auroc.dir/utils/parse_options.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/proxy_auroc.dir/utils/parse_options.cpp.o -c /Users/meteakgun/CLionProjects/CECILIA/utils/parse_options.cpp

CMakeFiles/proxy_auroc.dir/utils/parse_options.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/proxy_auroc.dir/utils/parse_options.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/meteakgun/CLionProjects/CECILIA/utils/parse_options.cpp > CMakeFiles/proxy_auroc.dir/utils/parse_options.cpp.i

CMakeFiles/proxy_auroc.dir/utils/parse_options.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/proxy_auroc.dir/utils/parse_options.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/meteakgun/CLionProjects/CECILIA/utils/parse_options.cpp -o CMakeFiles/proxy_auroc.dir/utils/parse_options.cpp.s

# Object files for target proxy_auroc
proxy_auroc_OBJECTS = \
"CMakeFiles/proxy_auroc.dir/apps/auroc/proxy.cpp.o" \
"CMakeFiles/proxy_auroc.dir/utils/parse_options.cpp.o"

# External object files for target proxy_auroc
proxy_auroc_EXTERNAL_OBJECTS =

proxy_auroc: CMakeFiles/proxy_auroc.dir/apps/auroc/proxy.cpp.o
proxy_auroc: CMakeFiles/proxy_auroc.dir/utils/parse_options.cpp.o
proxy_auroc: CMakeFiles/proxy_auroc.dir/build.make
proxy_auroc: CMakeFiles/proxy_auroc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/meteakgun/CLionProjects/CECILIA/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable proxy_auroc"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/proxy_auroc.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/proxy_auroc.dir/build: proxy_auroc

.PHONY : CMakeFiles/proxy_auroc.dir/build

CMakeFiles/proxy_auroc.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/proxy_auroc.dir/cmake_clean.cmake
.PHONY : CMakeFiles/proxy_auroc.dir/clean

CMakeFiles/proxy_auroc.dir/depend:
	cd /Users/meteakgun/CLionProjects/CECILIA/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/meteakgun/CLionProjects/CECILIA /Users/meteakgun/CLionProjects/CECILIA /Users/meteakgun/CLionProjects/CECILIA/cmake-build-debug /Users/meteakgun/CLionProjects/CECILIA/cmake-build-debug /Users/meteakgun/CLionProjects/CECILIA/cmake-build-debug/CMakeFiles/proxy_auroc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/proxy_auroc.dir/depend

