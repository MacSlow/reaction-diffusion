cmake_minimum_required (VERSION 2.8.12)
project (reaction-diffusion)

enable_language (ASM_NASM)
add_library (AsmLib src/laplace.nasm)

include (FindPkgConfig)
pkg_check_modules (SDL2 sdl2 REQUIRED)
include_directories (include ${SDL2_INCLUDE_DIRS})

add_compile_options(-std=c++14 -Werror -Wall -pedantic -Ofast -march=native)

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

add_custom_target (perf
	DEPENDS reaction-diffusion
	COMMAND perf stat ./reaction-diffusion
)

add_custom_target (report
	DEPENDS reaction-diffusion
	COMMAND perf record -g ./reaction-diffusion
	COMMAND perf report -g 'graph,0.5,caller' --sort comm,dso,sym
)
