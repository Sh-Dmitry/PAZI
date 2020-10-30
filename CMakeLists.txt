cmake_minimum_required(VERSION 3.5)

project(lab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_package(PkgConfig REQUIRED)
pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmpxx)

add_executable(lab main.cpp )


target_link_libraries(lab PRIVATE PkgConfig::gmp)
