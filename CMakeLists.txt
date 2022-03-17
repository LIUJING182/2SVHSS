cmake_minimum_required(VERSION 3.17)
project(ISC2019_Arb C)

set(CMAKE_C_STANDARD 99)

add_executable(ISC2019_Arb main.c PKE.h GaussRand.h VHSS.h GenData.h)

target_link_libraries(ISC2019_Arb libgmp.so  /usr/local/lib/libflint.so m.so)