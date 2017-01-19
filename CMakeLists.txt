################################################################################
#23456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
#
#  reaction-diffusion - 2D reaction-diffusion simulation pitting compiler-tuned
#                       vs hand-tuned x86/SIMD-assembly code
#
#  Copyright 2016-2017 Mirco Müller
# 
#  Author(s):
#    Mirco "MacSlow" Müller <macslow@gmail.com>
# 
#  This program is free software: you can redistribute it and/or modify it
#  under the terms of the GNU General Public License version 3, as published
#  by the Free Software Foundation.
# 
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranties of
#  MERCHANTABILITY, SATISFACTORY QUALITY, or FITNESS FOR A PARTICULAR
#  PURPOSE.  See the GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License along
#  with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
################################################################################

cmake_minimum_required (VERSION 2.8.12)
project (reaction-diffusion)

enable_language (ASM_NASM)
add_library (AsmLib src/laplace.nasm)

include (FindPkgConfig)
pkg_check_modules (SDL2 sdl2 REQUIRED)
include_directories (include ${SDL2_INCLUDE_DIRS})

add_compile_options(-std=c++14 -Werror -Wall -pedantic -fno-omit-frame-pointer -Ofast -march=native)

add_library (DisplayLib src/display.cpp)
add_library (BufferLib src/buffer.cpp)

add_executable (reaction-diffusion src/main.cpp)
target_link_libraries (reaction-diffusion
	DisplayLib
	BufferLib
	AsmLib
	${SDL2_LDFLAGS}
	-pthread
	-s
)

add_custom_target (valgrind
	DEPENDS reaction-diffusion
	COMMAND valgrind --track-origins=yes --show-leak-kinds=all  --leak-check=full -v ./reaction-diffusion
)

add_custom_target (perf
	DEPENDS reaction-diffusion
	COMMAND perf stat ./reaction-diffusion
)

add_custom_target (report
	DEPENDS reaction-diffusion
	COMMAND perf record -g ./reaction-diffusion
	COMMAND perf report -g 'graph,0.5,caller' --sort comm,dso,sym
)
