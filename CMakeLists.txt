cmake_minimum_required(VERSION 3.14)

project(svdc
    VERSION 1.0.0
    DESCRIPTION "Image compression with a truncated SVD written in pure C"
    LANGUAGES C)

option(SVDC_BUILD_TESTS "Build the unit tests" ON)
option(SVDC_OPENMP      "Compress color channels in parallel with OpenMP (if available)" ON)
option(SVDC_WERROR      "Treat compiler warnings as errors" OFF)
option(SVDC_SANITIZE    "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(SVDC_NATIVE      "Optimize for the build machine (-march=native)" OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------------------------------------------------------------
# warning flags for the project's own code (not for vendored stb)

if(MSVC)
    set(SVDC_WARNINGS /W4 /permissive-)
    if(SVDC_WERROR)
        list(APPEND SVDC_WARNINGS /WX)
    endif()
else()
    set(SVDC_WARNINGS
        -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wcast-qual
        -Wstrict-prototypes -Wmissing-prototypes -Wdouble-promotion -Wundef
        -Wformat=2)
    if(SVDC_WERROR)
        list(APPEND SVDC_WARNINGS -Werror)
    endif()
endif()

function(svdc_target_setup target)
    target_compile_options(${target} PRIVATE ${SVDC_WARNINGS})
    if(MSVC)
        # fopen and friends are standard C; silence the "use fopen_s" advice
        target_compile_definitions(${target} PRIVATE _CRT_SECURE_NO_WARNINGS)
    endif()
endfunction()

if(SVDC_SANITIZE AND NOT MSVC)
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all)
    add_link_options(-fsanitize=address,undefined)
endif()

if(SVDC_NATIVE AND NOT MSVC)
    add_compile_options(-march=native)
endif()

# libm is a separate library on most Unix systems
include(CheckLibraryExists)
check_library_exists(m sqrt "" SVDC_HAVE_LIBM)

# ---------------------------------------------------------------------------
# core library

add_library(svdcore STATIC
    src/matrix.c
    src/svd.c)
target_include_directories(svdcore PUBLIC src)
svdc_target_setup(svdcore)
if(SVDC_HAVE_LIBM)
    target_link_libraries(svdcore PUBLIC m)
endif()

# ---------------------------------------------------------------------------
# tests

if(SVDC_BUILD_TESTS)
    enable_testing()
    foreach(name test_svd)
        add_executable(${name} tests/${name}.c)
        target_link_libraries(${name} PRIVATE svdcore)
        svdc_target_setup(${name})
        add_test(NAME ${name} COMMAND ${name})
    endforeach()

    # matrix in, SVD out: used by tests/numpy_check.py
    add_executable(svd_tool tests/svd_tool.c)
    target_link_libraries(svd_tool PRIVATE svdcore)
    svdc_target_setup(svd_tool)

    # optional cross-check against numpy, registered only when numpy is importable
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy"
                        RESULT_VARIABLE SVDC_NUMPY_MISSING OUTPUT_QUIET ERROR_QUIET)
        if(SVDC_NUMPY_MISSING EQUAL 0)
            add_test(NAME numpy_check
                     COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/numpy_check.py
                             --tool $<TARGET_FILE:svd_tool>)
        endif()
    endif()
endif()
