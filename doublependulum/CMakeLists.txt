cmake_minimum_required(VERSION 3.17)
project(doublependulum2)

set(CMAKE_CXX_STANDARD 17)

add_executable(doublependulum2 main.cpp)

target_compile_definitions(${PROJECT_NAME} PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:_UNICODE UNICODE>)