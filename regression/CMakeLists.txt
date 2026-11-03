file(STRINGS "cases.txt" TEST_CASES)

include(esma_regression_tests)
esma_add_regression_tests(
  NAME fv3
  DATA_PATH FVdycoreCubed_GridComp
  TEST_CASES ${TEST_CASES}
)
