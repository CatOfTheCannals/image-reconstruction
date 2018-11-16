cmake_minimum_required(VERSION 3.5.1)
project(image-reconstruction CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -o3 -std=c++11 -Wall -Wno-long-long -pedantic -Wno-variadic-macros")

set(SOURCE_FILES
        src/Matrix.cpp
		src/Svd.cpp
	)
add_library(image-reconstruction-base ${SOURCE_FILES})

add_executable(main main.cpp)
target_link_libraries(main image-reconstruction-base)

add_subdirectory(test)

# Cosas importantes:
#   1) Google test genera su propio main, dejar main.cpp por afuera de SOURCE_FILES para evitar problemas
#   2) Al compilar "main", agregar el archivo main.cpp para que lo ejecute
