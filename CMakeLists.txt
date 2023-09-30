cmake_minimum_required(VERSION 3.26)
set(PROJECT_NAME gl_test)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 23)

add_compile_options(-march=native)
#add_compile_options(-mavx)
#add_compile_options(-mavx2)
#add_compile_options(-msse)

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

set(ASSIMP_BUILD_ALL_IMPORTERS_BY_DEFAULT OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_ALL_EXPORTERS_BY_DEFAULT OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_OBJ_IMPORTER ON CACHE BOOL "" FORCE)

add_subdirectory(libs/glfw-3.3.8)
add_subdirectory(libs/glm)
add_subdirectory(libs/assimp-master)

include_directories(libs/glad-2.0/include)
include_directories(libs/glm)
include_directories(libs/stbi/include)

add_executable(
    ${PROJECT_NAME} main.cpp
    libs/glad-2.0/src/glad.c
    libs/stbi/src/stb_image.cpp
        glEnv/linalg.hpp
        glEnv/shader.hpp
        glEnv/camera.hpp
        glEnv/primitives.hpp
        glEnv/lighting.hpp
        glEnv/geometry.hpp
        glEnv/utils.hpp
)

target_link_libraries(${PROJECT_NAME} glfw assimp)
