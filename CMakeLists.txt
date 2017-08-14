#
# @BEGIN LICENSE
#
# resp2 by Robert T. McGibbon, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(resp2 CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(psi4 1.1 REQUIRED)
find_package(NLOPT REQUIRED)
find_package(Boost REQUIRED)

add_psi4_plugin(resp2 src/dotsphere.cc
                      src/potential.cc
                      src/resp2.cc
                      src/respfit.cc
                      src/vdwsurface.cc)
target_include_directories(resp2 PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                 PRIVATE ${NLOPT_INCLUDE_DIRS}
                                 PRIVATE ${Boost_INCLUDE_DIRS})

target_link_libraries(resp2 PRIVATE "${NLOPT_LIBRARIES}")

# may be needed
