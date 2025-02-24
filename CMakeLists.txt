cmake_minimum_required(VERSION 3.10)
project(VideoTaps)

cmake_policy(SET CMP0091 NEW)
message(STATUS "Policy CMP0091 set NEW")
if(EXISTS "${CMAKE_BINARY_DIR}/generators/conan_toolchain.cmake")
    include("${CMAKE_BINARY_DIR}/generators/conan_toolchain.cmake")
endif()
add_compile_options(-w)

# Detectar el sistema operativo
if(WIN32)
    message(STATUS "Configurando para Windows")
    # Configuraciones específicas de Windows
    # Definir la macro _CRT_SECURE_NO_WARNINGS
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
elseif(UNIX)
    if(NOT CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "${CMAKE_BINARY_DIR}/Release/generators/conan_toolchain.cmake" CACHE STRING "Toolchain file")
    endif()
    message(STATUS "Configurando para Unix/Linux")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Release" FORCE)

    # Configuraciones específicas de Unix/Linux
endif()

# Define la ruta al archivo toolchain en el proyecto raíz
set(SENSIA_TOOLCHAIN_FILE "${CMAKE_BINARY_DIR}/generators/conan_toolchain.cmake" CACHE STRING "Toolchain file")

# Usando pkg_check_modules para buscar Aravis

# Habilitar C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Buscar SENSIANuc.h, .cpp
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

message("SRC_DIR: ${SRC_DIR}")


set(SOURCES
    ${SRC_DIR}/main.cpp
)

# Añadir rutas de Conan para encontrar paquetes
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}")

# Buscar Aravis
find_path(aravis_INCLUDE_DIR
    NAMES arv.h
    PATHS ${aravis_PKGCONF_INCLUDE_DIRS}
)
# Buscar OpenCV
find_path(OpenCV_INCLUDE_DIR
    NAMES opencv2/opencv.hpp
    PATHS ${OpenCV_PKGCONF_INCLUDE_DIRS}
)

# Buscar doctest
find_path(doctest_INCLUDE_DIR
    NAMES doctest.h
    PATHS ${doctest_PKGCONF_INCLUDE_DIRS}
)

# message("Conan Include Dirs: ${CONAN_INCLUDE_DIRS}")

find_package(aravis REQUIRED)
find_package(OpenCV REQUIRED)
find_package(doctest REQUIRED)


# Condicional para compilar bibliotecas
if(DEFINED LibRef)
    message(STATUS "Compiling with LibRef: ${LibRef}")

    if(LibRef STREQUAL "libAravis")
        add_subdirectory(Libs/libaravis)
        target_compile_definitions(${PROJECT_NAME} PUBLIC USE_PLUGIN)
        target_link_libraries(${PROJECT_NAME} PUBLIC libAravis)

    endif()
endif()

include_directories(${FOLDER_DIR})
# add_executable(${PROJECT_NAME} src/main.cpp ../SensiaNuc.cpp ../pixel_manager.cpp)
add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ${CONAN_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PUBLIC ${CONAN_LIBS} ${OpenCV_LIBS} aravis::aravis)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)



