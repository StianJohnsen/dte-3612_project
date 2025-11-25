cmake_minimum_required(VERSION 3.19)
project(fem_implementation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



add_executable(fem_implementation
main.cpp
FemObject.h
FemObject.cpp
)


set(ttl_DIR "C:/Program Files/SINTEF/ttl/lib/cmake/ttl")
find_package(ttl CONFIG REQUIRED)
# find_package(ttl CONFIG REQUIRED)

# include_directories("C:/dev/skole/dte-3612/TTL/include")

find_package(Eigen3 CONFIG REQUIRED)

target_link_libraries(fem_implementation PRIVATE
    ttl
    Eigen3::Eigen
)
