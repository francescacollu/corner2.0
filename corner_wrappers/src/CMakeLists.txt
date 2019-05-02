set(CMAKE_CXX_FLAGS "-std=c++11")

find_package(Armadillo REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(corner_wrappers armadilloVector.cpp)

target_include_directories(corner_wrappers PRIVATE ../include/)
target_include_directories(corner_wrappers PRIVATE ${Eigen3_INCLUDE_DIR})
target_include_directories(corner_wrappers PRIVATE ${ARMADILLO_INCLUDE_DIRS})

target_link_libraries(corner_wrappers ${Armadillo_LIBRARIES})