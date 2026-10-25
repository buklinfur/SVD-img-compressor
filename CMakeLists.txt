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
    set(SVDC_WARNINGS /W4)
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
# vendored stb: compiled in its own object library without our warning flags

add_library(svdc_stb OBJECT src/stb_impl.c)
target_include_directories(svdc_stb SYSTEM PUBLIC third_party/stb)
if(MSVC)
    target_compile_options(svdc_stb PRIVATE /w)
    target_compile_definitions(svdc_stb PRIVATE _CRT_SECURE_NO_WARNINGS)
else()
    target_compile_options(svdc_stb PRIVATE -w)
    if(SVDC_SANITIZE)
        # stb_image_write's JPEG encoder has a signed left-shift overflow (upstream UB,
        # harmless on every real compiler); keep all other checks for stb
        target_compile_options(svdc_stb PRIVATE -fno-sanitize=shift)
    endif()
endif()
set_target_properties(svdc_stb PROPERTIES POSITION_INDEPENDENT_CODE ON)

# ---------------------------------------------------------------------------
# core library

add_library(svdcore STATIC
    src/matrix.c
    src/svd.c
    src/image.c
    src/codec.c
    src/metrics.c
    $<TARGET_OBJECTS:svdc_stb>)
target_include_directories(svdcore PUBLIC src)
target_include_directories(svdcore SYSTEM PUBLIC third_party/stb)
svdc_target_setup(svdcore)
if(SVDC_HAVE_LIBM)
    target_link_libraries(svdcore PUBLIC m)
endif()

if(SVDC_OPENMP AND MSVC)
    # MSVC's /openmp (OpenMP 2.0) cannot be combined with /std:c11
    message(STATUS "svdc: OpenMP is not used with MSVC in C11 mode, channels are processed sequentially")
elseif(SVDC_OPENMP)
    find_package(OpenMP COMPONENTS C)
    if(OpenMP_C_FOUND)
        target_link_libraries(svdcore PUBLIC OpenMP::OpenMP_C)
        message(STATUS "svdc: OpenMP enabled")
    else()
        message(STATUS "svdc: OpenMP not found, channels are processed sequentially")
    endif()
endif()

# ---------------------------------------------------------------------------
# command line tool

add_executable(svdc src/main.c)
target_link_libraries(svdc PRIVATE svdcore)
target_compile_definitions(svdc PRIVATE SVDC_VERSION_STRING="${PROJECT_VERSION}")
svdc_target_setup(svdc)

include(GNUInstallDirs)
install(TARGETS svdc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# ---------------------------------------------------------------------------
# tests

if(SVDC_BUILD_TESTS)
    enable_testing()
    foreach(name test_svd test_codec)
        add_executable(${name} tests/${name}.c)
        target_link_libraries(${name} PRIVATE svdcore)
        target_compile_definitions(${name} PRIVATE EXAMPLES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/examples")
        svdc_target_setup(${name})
        add_test(NAME ${name} COMMAND ${name})
    endforeach()

    # end-to-end tests of the command line tool
    set(EX ${CMAKE_CURRENT_SOURCE_DIR}/examples)
    set(OUT ${CMAKE_CURRENT_BINARY_DIR}/cli-test)
    file(MAKE_DIRECTORY ${OUT})
    add_test(NAME cli_compress_rgba COMMAND svdc compress ${EX}/test3.png ${OUT}/rgba.svdc -k 30)
    add_test(NAME cli_decompress_rgba COMMAND svdc decompress ${OUT}/rgba.svdc ${OUT}/rgba.png)
    add_test(NAME cli_compress_jpeg_energy COMMAND svdc compress ${EX}/img.jpg ${OUT}/jpeg.svdc --energy 0.99 --bits 16 -q)
    add_test(NAME cli_decompress_to_jpeg COMMAND svdc decompress ${OUT}/jpeg.svdc ${OUT}/jpeg.jpg -q)
    add_test(NAME cli_info COMMAND svdc info ${OUT}/rgba.svdc)
    add_test(NAME cli_jacobi COMMAND svdc compress ${EX}/test.png ${OUT}/small.svdc --method jacobi -k 100)
    add_test(NAME cli_bench COMMAND svdc bench ${EX}/test.png --ks 1,10,40 --csv ${OUT}/bench.csv)
    set_tests_properties(cli_compress_rgba PROPERTIES FIXTURES_SETUP rgba)
    set_tests_properties(cli_decompress_rgba cli_info PROPERTIES FIXTURES_REQUIRED rgba)
    set_tests_properties(cli_compress_jpeg_energy PROPERTIES FIXTURES_SETUP jpeg)
    set_tests_properties(cli_decompress_to_jpeg PROPERTIES FIXTURES_REQUIRED jpeg)

    # these must fail cleanly (non-zero exit, no crash)
    add_test(NAME cli_rejects_non_svdc COMMAND svdc decompress ${EX}/test.png ${OUT}/x.png)
    add_test(NAME cli_rejects_missing_file COMMAND svdc compress ${EX}/missing.png ${OUT}/x.svdc)
    add_test(NAME cli_rejects_bad_k COMMAND svdc compress ${EX}/test.png ${OUT}/x.svdc -k 0)
    add_test(NAME cli_rejects_k_and_energy COMMAND svdc compress ${EX}/test.png ${OUT}/x.svdc -k 5 --energy 0.9)
    set_tests_properties(cli_rejects_non_svdc cli_rejects_missing_file cli_rejects_bad_k
                         cli_rejects_k_and_energy PROPERTIES WILL_FAIL TRUE)

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
