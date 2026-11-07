cmake_minimum_required(VERSION 3.15)
project(hydrogen_orbital_explorer VERSION 1.0.0 LANGUAGES CXX)

option(BUILD_VISUALIZER "Build the macOS OpenGL application" ON)
option(BUILD_TOOLS "Build the headless point-cloud exporter" ON)
include(CTest)

add_library(quantum_core
    src/quantum.cpp
    src/sampling.cpp
    src/flow.cpp
)
target_include_directories(quantum_core PUBLIC src)
target_compile_features(quantum_core PUBLIC cxx_std_17)
set_target_properties(quantum_core PROPERTIES CXX_EXTENSIONS OFF)

function(enable_project_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endfunction()
enable_project_warnings(quantum_core)

if(BUILD_TOOLS)
    add_executable(export_cloud tools/export_cloud.cpp)
    target_link_libraries(export_cloud PRIVATE quantum_core)
    set_target_properties(export_cloud PROPERTIES CXX_EXTENSIONS OFF)
    enable_project_warnings(export_cloud)
endif()

if(BUILD_TESTING)
    add_executable(quantum_tests tests/quantum_tests.cpp)
    target_link_libraries(quantum_tests PRIVATE quantum_core)
    set_target_properties(quantum_tests PROPERTIES CXX_EXTENSIONS OFF)
    enable_project_warnings(quantum_tests)
    add_test(NAME quantum_tests COMMAND quantum_tests)
endif()

if(BUILD_VISUALIZER)
    if(NOT APPLE)
        message(FATAL_ERROR
            "The OpenGL visualizer currently supports macOS only. "
            "Use -DBUILD_VISUALIZER=OFF for the dependency-free core and tests.")
    endif()
    # Homebrew defaults; explicit CMAKE_PREFIX_PATH and GLM_INCLUDE_DIR also work.
    list(APPEND CMAKE_PREFIX_PATH /opt/homebrew /usr/local)
    find_package(glfw3 3.3 CONFIG REQUIRED)
    find_path(GLM_INCLUDE_DIR glm/glm.hpp)
    if(NOT GLM_INCLUDE_DIR)
        message(FATAL_ERROR "GLM headers not found. On macOS: brew install glfw glm")
    endif()
    find_package(OpenGL REQUIRED)
    find_package(Threads REQUIRED)

    add_executable(quantum_orbitals src/main.cpp src/renderer.cpp)
    target_include_directories(quantum_orbitals PRIVATE ${GLM_INCLUDE_DIR})
    target_link_libraries(quantum_orbitals PRIVATE
        quantum_core glfw OpenGL::GL Threads::Threads)
    target_compile_definitions(quantum_orbitals PRIVATE GL_SILENCE_DEPRECATION)
    set_target_properties(quantum_orbitals PROPERTIES CXX_EXTENSIONS OFF)
    enable_project_warnings(quantum_orbitals)
endif()
