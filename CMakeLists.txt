cmake_minimum_required(VERSION 2.8)
project(SimpleDGLAP CXX)


include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++17" CXX17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

link_libraries(stdc++fs)

if(NOT CXX17)
  message(FATAL_ERROR "Your compiler does not support C++17")
endif()

# default build type: Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# default install prefix: ~/.local
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local"
      CACHE PATH "Install path prefix, prepended onto install directories."
      FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

#string(APPEND CMAKE_CXX_FLAGS " -march=native")
# By default, optimize for the system's native architecture.  Disable via the
# NATIVE option.  In addition, detect if another architecture flag is already
# set and do not override it.
option(NATIVE "compile for native architecture" ON)
if(NATIVE)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # Intel compiler: search for -m, -x, -ax flags; if not found add -xHost.
    if(NOT CMAKE_CXX_FLAGS MATCHES "(^| )-(m|x|ax)[^ ]+")
      string(APPEND CMAKE_CXX_FLAGS " -xHost")
    endif()
  else()
    # All other compilers (gcc-like): search for -m flags;
    # if not found add -march=native.
    if(NOT CMAKE_CXX_FLAGS MATCHES "(^| )-m[^ ]+")
      string(APPEND CMAKE_CXX_FLAGS " -march=native")
    endif()
  endif()
endif()

#string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra")

option(MORE_WARNINGS "enable more compiler warnings" OFF)
if(MORE_WARNINGS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # adapted from http://stackoverflow.com/a/9862800
  string(APPEND CMAKE_CXX_FLAGS " -Wpedantic -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-conversion -Wsign-promo -Wstrict-overflow=2 -Wundef")
endif()

# disable silly clang warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  string(APPEND CMAKE_CXX_FLAGS " -Wno-missing-braces -Wno-c++17-narrowing")
endif()

string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Werror")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Og")
endif()


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  # Check for sanitizer support.
  # (see e.g. https://code.google.com/p/address-sanitizer/wiki/AddressSanitizer)
  # check_cxx_compiler_flag() doesn't work here because it doesn't set the test
  # flag at link time, so use check_cxx_source_compiles() instead.
  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_QUIET TRUE)

  # First try both AddressSanitizer and UndefinedBehaviorSanitizer.
  # Disable the vptr sanitizer because it throws errors with boost::any.
  set(CMAKE_REQUIRED_FLAGS "-fsanitize=address,undefined -fno-sanitize=vptr")
  check_cxx_source_compiles("int main() { return 0; }" AddrUndefSanitizer)

  # Lacking that, try AddressSanitizer only.
  if(NOT AddrUndefSanitizer)
    set(CMAKE_REQUIRED_FLAGS "-fsanitize=address")
    check_cxx_source_compiles("int main() { return 0; }" AddrSanitizer)
    if(NOT AddrSanitizer)
      # out of luck...
      set(CMAKE_REQUIRED_FLAGS "")
    endif()
  endif()
endif()

string(APPEND CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_REQUIRED_FLAGS} --coverage")

set(LIBRARY_NAME "lib${PROJECT_NAME}")
add_subdirectory(src)
include_directories(src)


foreach(App "light" "charm" "bottom")
	add_executable(${App} ./TestMains/${App}.cpp)
	target_link_libraries(${App} ${LIBRARY_NAME})
	install(TARGETS ${App} DESTINATION bin)
endforeach()

