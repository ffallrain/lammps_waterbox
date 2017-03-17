cmake_minimum_required(VERSION 3.4)
project(pro1)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++11 -nostdinc")
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)

include_directories("/usr/include"
        )

add_executable("pro1" main.cpp)


#add_custom_target(project COMMAND make serial -C /home/zengping/work/code/waterbox/lammps-master//src
#        )
add_custom_target(project COMMAND make serial -C ${PROJECT_SOURCE_DIR}/src
        )