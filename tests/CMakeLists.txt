# Minima versione di CMake richiesta per buildare il tuo progetto
# Serve almeno la 2.8.12 per usare i comandi target_link_libraries ecc...
cmake_minimum_required(VERSION 2.8.12)

# Digli di usare lo standard c++11
set(CMAKE_CXX_STANDARD 11)

find_package(Armadillo REQUIRED)
message("ARMADILLO LIBRARY PATH: " ${ARMADILLO_LIBRARIES})

# Setta le flag per compilare root da dare al compilatore
set(CMAKE_CXX_FLAGS ${ROOT_CXX_FLAGS})
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Qua gli dici come si chiama il tuo eseguibile
# e quali sorgenti deve compilare
# Dato che EXEC_NAME � uguale a corner2.0
# l'eseguibile si chiamer� cos�
add_executable(testSparse testSparse.cpp)
add_executable(testSort testSort.cpp
    ../src/eigenV.cpp
    ../src/utility.cpp)


# Aggiungi la cartella degli header di Root e Armadillo alla
# lista di cartelle in cui guardare quando dal codice tu
# scrivi #include ...
# In questo modo quando scrivi #include <armadillo> lui guarda
# anche nella cartella di Armadillo
target_include_directories(testSparse PUBLIC ${ROOT_INCLUDE_DIRS})
target_include_directories(testSparse PUBLIC ${ARMADILLO_INCLUDE_DIRS})

target_include_directories(testSort PUBLIC ../include)
target_include_directories(testSort PUBLIC ${ARMADILLO_INCLUDE_DIRS})


message("Armadillo_include: ${ARMADILLO_INCLUDE_DIRS}")

# Linka le librerie al tuo eseguibile
target_link_libraries(testSparse ${ROOT_ALL_LIBRARIES} ${ARMADILLO_LIBRARIES} ${ARPACK_LIBRARY})
target_link_libraries(testSort ${ARMADILLO_LIBRARIES})