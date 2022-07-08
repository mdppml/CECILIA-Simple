cmake_minimum_required(VERSION 3.9)
project(CECILIA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_DEBUG "-O3")
set(CMAKE_CXX_FLAGS_RELEASE  "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(BENCHMARK_ENABLE_TESTING OFF) # required to disable google benchmark's google testing dependency

add_executable(helper apps/helper.cpp)

#add_executable(proxy_auroc apps/auroc/proxy.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h utils/llib.h ${ENC_FILES})

#add_executable(proxy_auroc_tie apps/auroctie/proxy.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h ${ENC_FILES})

#add_executable(proxy_aupr apps/aupr/proxy.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h utils/llib.h ${ENC_FILES} core/auc.h)

add_executable(helper_cnn apps/cnn/helper.cpp core/Party.h core/core.h core/cnn.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h utils/llib.h ${ENC_FILES} apps/cnn/model_parser.h )
add_executable(proxy_cnn apps/cnn/proxy.cpp core/Party.h core/core.h core/cnn.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h utils/llib.h ${ENC_FILES} apps/cnn/model_parser.h )

add_executable(proxy_test apps/test/proxy.cpp core/Party.h core/core.h core/cnn.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h utils/llib.h ${ENC_FILES})

add_executable(proxy_benchmark apps/benchmark/proxy.cpp)

add_executable(proxy_rkn apps/rkn/proxy.cpp core/Party.h core/rkn.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h ${ENC_FILES})

# include dependencies using git submodule:
add_subdirectory(dependencies/eigen)
add_subdirectory(dependencies/benchmark)

# link libraries:
find_package(Threads REQUIRED)

target_link_libraries(helper PRIVATE Threads::Threads)

target_link_libraries(proxy_test PRIVATE Threads::Threads)

target_link_libraries(helper_cnn PRIVATE Threads::Threads)
target_link_libraries(proxy_cnn PRIVATE Threads::Threads)

target_link_libraries(proxy_benchmark PRIVATE Threads::Threads benchmark::benchmark)

target_link_libraries(proxy_rkn PRIVATE Threads::Threads Eigen3::Eigen)