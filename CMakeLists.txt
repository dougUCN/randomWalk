cmake_minimum_required(VERSION 2.8)
project(randomWalk_t)

set(CMAKE_CXX_STANDARD 11)
include_directories("include")

set(Boost_REALPATH ON)
find_package(Boost REQUIRED COMPONENTS program_options)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
else()
    message(SEND_ERROR "Boost not found!")
endif()

find_package(HDF5 REQUIRED COMPONENTS CXX HL)
if (HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
else()
    message(SEND_ERROR "HDF5 not found! Try: sudo apt-get install -y libhdf5-dev")
endif()

# List of executables
add_executable( randomWalk_t.x src/main.cpp src/particle1d.cpp src/mc.cpp src/output.cpp)
target_link_libraries( randomWalk_t.x
                        ${Boost_LIBRARIES}
                        ${HDF5_CXX_LIBRARIES}
                        ${HDF5_HL_LIBRARIES}
                        ${HDF5_LIBRARIES})
