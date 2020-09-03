set(BINARY run_tests)

set(TEST_SRCS
  basic.cpp )

ADD_EXECUTABLE(${BINARY} ${TEST_SRCS})

target_link_libraries(${BINARY} gtest gtest_main)

add_test(AllTests ${BINARY})

