set(TEST_BINARY run_tests)

set(TEST_FILES
  all_tests.cpp )

ADD_EXECUTABLE(${TEST_BINARY} ${TEST_FILES})

target_link_libraries(${TEST_BINARY} PUBLIC libSimPLEX gtest gtest_main)

add_test(AllTests ${TEST_BINARY})

