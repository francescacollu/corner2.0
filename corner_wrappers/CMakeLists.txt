cmake_minimum_required(VERSION 2.6)
project(corner_wrappers)

add_subdirectory(tests)
add_subdirectory(src)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake")