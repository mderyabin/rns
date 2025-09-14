enable_testing()

include(GoogleTest)

file (GLOB TESTS_SRC_FILES CONFIGURE_DEPENDS src/*.cpp)
foreach (app ${TESTS_SRC_FILES})
    get_filename_component ( exe ${app} NAME_WE )
    add_executable ( ${exe} ${app} )
    set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DERICTORY}/tests)
    target_link_libraries(${exe} PRIVATE modular)
    target_link_libraries(${exe} PUBLIC gmp)
    target_link_libraries(${exe} PUBLIC m)
    # target_link_libraries(${exe} PRIVATE GTest::gtest)
    target_link_libraries(${exe} PRIVATE GTest::gtest_main)

    gtest_discover_tests(${exe})
endforeach()
