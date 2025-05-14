cmake_minimum_required(VERSION 3.10)

# Projektname und Version
project(lorentz_gas_simulator VERSION 1.0)

# Setze den Compiler auf mpic++
set(CMAKE_CXX_COMPILER mpic++)

# C++-Standard festlegen
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include-Pfade hinzufügen
include_directories(${PROJECT_SOURCE_DIR}/include)

# Quellen und Header-Dateien sammeln
file(GLOB_RECURSE SOURCES 
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
    "${PROJECT_SOURCE_DIR}/include/*.h"
)

# Ausführbare Datei erstellen
add_executable( lorentz_mpi_0 src/main_mpi_0.cpp src/globals.cpp )
add_executable( lorentz_mpi_B src/main_mpi_B.cpp src/globals.cpp )

# Bibliotheken verlinken (falls vorhanden)
# target_link_libraries(${PROJECT_NAME} PRIVATE your_library)

# Tests hinzufügen (optional, falls eine Testdatei im tests/ Ordner vorhanden ist)
file(GLOB_RECURSE TEST_SOURCES "${PROJECT_SOURCE_DIR}/tests/*.cpp")
if(TEST_SOURCES)
    enable_testing()
    add_executable(run_tests ${TEST_SOURCES})
    add_test(NAME run_tests COMMAND run_tests)
endif()

# Installationseinrichtung (optional)
install(TARGETS lorentz_mpi_0 lorentz_mpi_B DESTINATION bin)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)
