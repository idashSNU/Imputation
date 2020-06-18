set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
cmake_minimum_required(VERSION 3.10)

include_directories(../../ModHEaaN/heaan/include/)



link_directories(../../ModHEaaN/lib/)

set (CMAKE_CXX_STANDARD 14)

add_compile_options(-Ofast)

add_compile_options(-fopenmp)

add_compile_options(-pthread)


ADD_EXECUTABLE( enc_impute
                enc_impute.cpp
                ../src/TimeUtils.cpp
                ../src/PrintUtils.cpp
				../src/Evaluator.cpp
				../src/FileReader.cpp
                ../src/SNP_Parameters.cpp)

TARGET_LINK_LIBRARIES( enc_impute
                       HEAAN pthread crypto gomp)
