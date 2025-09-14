file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cpp)

add_library(modular STATIC ${SRC_FILES})

target_include_directories(modular PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DERICTORY ${CMAKE_BINARY_DIR}/lib)

target_link_libraries(modular PUBLIC gmp)
target_link_libraries(modular PUBLIC m)
target_link_libraries(modular PUBLIC ntl)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)