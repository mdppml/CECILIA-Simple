cmake_minimum_required(VERSION 3.9)
project(CECILIA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_DEBUG "-O3")
set(CMAKE_CXX_FLAGS_RELEASE  "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


add_executable(helper_auroc apps/auroc/helper.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h ${ENC_FILES})
add_executable(proxy_auroc apps/auroc/proxy.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h apps/auroc/llib.h ${ENC_FILES})

add_executable(helper_auroc_tie apps/auroctie/helper.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h ${ENC_FILES})
add_executable(proxy_auroc_tie apps/auroctie/proxy.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h apps/auroctie/llib.h ${ENC_FILES})

add_executable(helper_aupr apps/aupr/helper.cpp core/Party.h  utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h ${ENC_FILES} core/auc.h)
add_executable(proxy_aupr apps/aupr/proxy.cpp core/Party.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h apps/aupr/llib.h ${ENC_FILES} core/auc.h)

add_executable(helper_test apps/test/helper.cpp core/Party.h core/core.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h ${ENC_FILES})
add_executable(proxy_test apps/test/proxy.cpp core/Party.h core/core.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h apps/test/llib.h ${ENC_FILES})

add_executable(helper_cnn apps/cnn/helper.cpp core/Party.h core/cnn.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h ${ENC_FILES})
add_executable(proxy_cnn apps/cnn/proxy.cpp core/Party.h core/cnn.h utils/constant.h utils/parse_options.cpp utils/parse_options.h utils/connection.h utils/flib.h apps/cnn/llib.h ${ENC_FILES})