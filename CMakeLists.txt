cmake_minimum_required(VERSION 3.16)
project(TurbDrive LANGUAGES CXX)

# -------------------------------------------------
# C++ standard
# -------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -------------------------------------------------
# Find dependencies
# -------------------------------------------------
find_package(FFTW3 REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C HL)

# -------------------------------------------------
# Executable
# -------------------------------------------------
add_executable(turbdrive
    main.cpp
    TurbDrive.cpp
    TurbDrive.h
)

# -------------------------------------------------
# Include directories
# -------------------------------------------------
target_include_directories(turbdrive PRIVATE
    ${FFTW3_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

# -------------------------------------------------
# Link libraries
# Use only absolute paths if needed (Homebrew fallback)
# -------------------------------------------------
# Default Homebrew paths (Apple Silicon)
if(NOT DEFINED FFTW3_LIB_PATH)
    set(FFTW3_LIB_PATH "/opt/homebrew/lib/libfftw3.dylib")
endif()
if(NOT DEFINED HDF5_LIB_PATH)
    set(HDF5_LIB_PATH "/opt/homebrew/lib/libhdf5.dylib")
endif()

target_link_libraries(turbdrive PRIVATE
    ${FFTW3_LIB_PATH}
    ${HDF5_LIB_PATH}
)
