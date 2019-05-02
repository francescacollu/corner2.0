set(CMAKE_CXX_FLAGS "-std=c++11")

find_package(Armadillo REQUIRED)
find_package(Eigen3 REQUIRED)

add_executable(testArmadilloVector testArmadilloVector.cpp)
target_include_directories(testArmadilloVector PRIVATE ../include ${ARMADILLO_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS})
message("Armadilo include dirs: ${ARMADILLO_INCLUDE_DIRS}")
target_link_libraries(testArmadilloVector corner_wrappers ${ARMADILLO_LIBRARIES})
