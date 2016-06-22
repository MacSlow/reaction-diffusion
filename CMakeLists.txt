cmake_minimum_required (VERSION 2.8.12)
project (reaction-diffusion)

add_compile_options(-std=c++14 -Werror -ggdb -Wall -pedantic -O3)
#add_compile_options(-std=c++14 -Werror -Wall -pedantic -O3 -ggdb)

include (FindPkgConfig)
pkg_check_modules (SDL2 sdl2 REQUIRED)

include_directories (include ${SDL2_INCLUDE_DIRS})

add_library (DisplayLib src/display.cpp)
add_library (BufferLib src/buffer.cpp)

set (asm_obj_files laplace.o)
add_library (AsmLib STATIC EXCLUDE_FROM_ALL ${asm_obj_files})
add_custom_command (OUTPUT laplace.o
	COMMAND nasm -f elf64 ../src/laplace.s -o laplace.o)
set_target_properties (AsmLib
	PROPERTIES
	LINKER_LANGUAGE C
	ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable (reaction-diffusion src/main.cpp)
target_link_libraries (reaction-diffusion
	DisplayLib
	BufferLib
	AsmLib
	${SDL2_LDFLAGS}
	-pthread
	-s
)
