# Nome del progetto
project(Tesi)


# Setta una variabile di nome EXEC_NAME uguale a StructuredCorner
set(EXEC_NAME corner2.0)
SET(ROOT_PATH "" CACHE STRING "path of ROOT base folder")
list(APPEND CMAKE_PREFIX_PATH ${ROOT_PATH}/cmake)

SET(ARMADILLO_PATH "" CACHE STRING "path of Armadillo base folder")
list(APPEND CMAKE_PREFIX_PATH ${ARMADILLO_PATH})

SET(ARPACK_LIBRARY "" CACHE STRING "path of Arpack library")

# Minima versione di CMake richiesta per buildare il tuo progetto
# Serve almeno la 2.8.12 per usare i comandi target_link_libraries ecc...
cmake_minimum_required(VERSION 2.8.12)

# Digli di usare lo standard c++11
set(CMAKE_CXX_STANDARD 11)


# Trova ROOT e ARMADILLO e definisce delle variabili tipo ROOT_INCLUDE_DIRS ecc...
# che ti servono per linkarle al tuo eseguibile. Se non lo trova ti d� errore
# perch� gli ho detto che sono REQUIRED
find_package(ROOT REQUIRED COMPONENTS Graf Gpad)

set(ROOT_ALL_LIBRARIES ${ROOT_LIBRARIES})

find_package(Armadillo REQUIRED)
message("ARMADILLO LIBRARY PATH: " ${ARMADILLO_LIBRARIES})

# Setta le flag per compilare root da dare al compilatore
set(CMAKE_CXX_FLAGS ${ROOT_CXX_FLAGS})
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Questa � la lista delle sorgenti da compilare per creare
# il tuo eseguibile
# Se ti d� fastidio scrivere ogni volta src
# c'� un modo per evitarlo ma non volevo mettere troppa roba
set(SOURCES
	src/testSystem.cpp
	src/Block.cpp
	src/HMatrix.cpp
	src/Site.cpp
	src/System.cpp
	src/utility.cpp
	src/eigenV.cpp
	# aggiungi qui tutte le sorgenti...	
)

# Qua gli dici come si chiama il tuo eseguibile
# e quali sorgenti deve compilare
# Dato che EXEC_NAME � uguale a corner2.0
# l'eseguibile si chiamer� cos�
add_executable(${EXEC_NAME} ${SOURCES})
add_executable(quantumTrajectories src/quantumTrajectories.cpp)
add_executable(plotter src/plotter.cpp src/utility.cpp)
add_executable(qTraj4sites src/qTraj4sites.cpp)


# Aggiungi la cartella degli header di Root e Armadillo alla
# lista di cartelle in cui guardare quando dal codice tu
# scrivi #include ...
# In questo modo quando scrivi #include <armadillo> lui guarda
# anche nella cartella di Armadillo
target_include_directories(corner2.0 PUBLIC ${ROOT_INCLUDE_DIRS})
target_include_directories(corner2.0 PUBLIC ${ARMADILLO_INCLUDE_DIRS})

target_include_directories(quantumTrajectories PUBLIC ${ROOT_INCLUDE_DIRS})
target_include_directories(quantumTrajectories PUBLIC ${ARMADILLO_INCLUDE_DIRS})

target_include_directories(qTraj4sites PUBLIC ${ROOT_INCLUDE_DIRS})
target_include_directories(qTraj4sites PUBLIC ${ARMADILLO_INCLUDE_DIRS})

target_include_directories(plotter PUBLIC ${ROOT_INCLUDE_DIRS})
target_include_directories(plotter PUBLIC ${ARMADILLO_INCLUDE_DIRS})

# Linka le librerie al tuo eseguibile
target_link_libraries(corner2.0 ${ROOT_ALL_LIBRARIES} ${ARMADILLO_LIBRARIES})
target_link_libraries(quantumTrajectories ${ROOT_ALL_LIBRARIES} ${ARMADILLO_LIBRARIES})
target_link_libraries(qTraj4sites ${ROOT_ALL_LIBRARIES} ${ARMADILLO_LIBRARIES})
target_link_libraries(plotter ${ROOT_ALL_LIBRARIES} ${ARMADILLO_LIBRARIES})

message("ROOT libraries path: ${ROOT_ALL_LIBRARIES}")

# Questa � un po' sottile. Per ora sappi solo che stiamo aggiungendo la cartella include
# alla lista delle cartelle in cui guarda il compilatore
# quando scrivi #include ..., come nel caso dell'altro comando usato prima
target_include_directories(${EXEC_NAME} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_include_directories(plotter PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

add_subdirectory(tests)
