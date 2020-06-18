# CMake instructions to make the static lib
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 14)

add_compile_options(-Ofast)

add_compile_options(-fPIC)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ../lib)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ../lib)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

include_directories(./include/)

link_directories(../lib)


FILE(GLOB srcs ./src/*.cpp)
FILE(GLOB basic_srcs ./src/basic/*.cpp)

ADD_LIBRARY(HEAAN SHARED ${basic_srcs} ${srcs})

target_link_libraries(HEAAN HEAAN pthread crypto)


install(DIRECTORY ./include DESTINATION .)
install(FILES ../lib/heaan.cpython-36m-x86_64-linux-gnu.so ../lib/libHEAAN.so DESTINATION lib)
