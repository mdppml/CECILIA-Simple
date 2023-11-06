cmake_minimum_required(VERSION 3.21)
project(CECILIA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE  "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Libraries:
find_package(Threads REQUIRED)


# Helper
add_executable(helper apps/helper.cpp)
target_link_libraries(helper PRIVATE Threads::Threads cryptopp::cryptopp)

# CNN
add_executable(helper_cnn apps/cnn/helper.cpp core/Party.h core/core.h core/cnn.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h utils/llib.h ${ENC_FILES} apps/cnn/model_parser.h )
add_executable(proxy_cnn apps/cnn/proxy.cpp core/Party.h core/core.h core/cnn.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h utils/llib.h ${ENC_FILES} apps/cnn/model_parser.h apps/cnn/mnist_reader.hpp apps/cnn/mnist_reader_common.hpp apps/cnn/mnist_utils.hpp)
target_link_libraries(helper_cnn PRIVATE Threads::Threads cryptopp::cryptopp)
target_link_libraries(proxy_cnn PRIVATE Threads::Threads cryptopp::cryptopp)

# CECILIA - ppRKN
add_executable(proxy_rkn apps/rkn/proxy.cpp core/Party.h core/rkn.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h ${ENC_FILES})
target_link_libraries(proxy_rkn PRIVATE Threads::Threads Eigen3::Eigen cryptopp::cryptopp)

# ppAURORA
add_executable(proxy_auroc apps/auroc/proxy.cpp core/Party.h utils/constant.h utils/auc_utils.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h utils/llib.h ${ENC_FILES})
target_link_libraries(proxy_auroc PRIVATE Threads::Threads cryptopp::cryptopp)
add_executable(proxy_auroc_tie apps/auroctie/proxy.cpp core/Party.h utils/constant.h utils/auc_utils.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h ${ENC_FILES})
target_link_libraries(proxy_auroc_tie PRIVATE Threads::Threads cryptopp::cryptopp)
add_executable(proxy_aupr apps/aupr/proxy.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h utils/llib.h ${ENC_FILES} core/auc.h)
target_link_libraries(proxy_aupr PRIVATE Threads::Threads cryptopp::cryptopp)

# Test
add_executable(proxy_test apps/test/mainTest.cpp core/Party.h core/core.h core/cnn.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h utils/llib.h ${ENC_FILES})
target_link_libraries(proxy_test PRIVATE Threads::Threads cryptopp::cryptopp)

# Unit Tests
add_executable(test apps/test/mainTest.cpp core/Party.h core/core.h core/cnn.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h utils/llib.h ${ENC_FILES})
target_link_libraries(test PRIVATE Threads::Threads cryptopp::cryptopp)

# Benchmark
add_executable(benchmark apps/benchmark/main.cpp)
target_link_libraries(benchmark PRIVATE Threads::Threads cryptopp::cryptopp)

# include dependencies using git submodule:
add_subdirectory(dependencies/eigen)
add_subdirectory(dependencies/cryptopp-cmake)
