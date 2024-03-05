cmake_minimum_required(VERSION 3.21)
project(CECILIA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_DEBUG "-O3")
set(CMAKE_CXX_FLAGS_RELEASE  "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Libraries:
find_package(Threads REQUIRED)

# Helper
add_executable(helper apps/helper.cpp)
target_link_libraries(helper PRIVATE Threads::Threads cryptopp::cryptopp)

# Demo
add_executable(demo apps/demo/proxy.cpp booleancore/core.h core/Party.h core/core.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h utils/llib.h ${ENC_FILES})
target_link_libraries(demo PRIVATE Threads::Threads cryptopp::cryptopp)

# Demo
add_executable(exercise apps/exercise/proxy.cpp booleancore/core.h core/Party.h core/core.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h utils/llib.h ${ENC_FILES})
target_link_libraries(exercise PRIVATE Threads::Threads cryptopp::cryptopp)


# include dependencies using git submodule:
add_subdirectory(dependencies/eigen)
add_subdirectory(dependencies/cryptopp-cmake)
