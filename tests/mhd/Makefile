all:
	@make tests/mhd_executables -C ../..

t: test
test:
	@make tests/mhd_tests -C ../..

.PHONY: results
r: results
results:
	@make tests/mhd_clean_results -C ../..

c: clean
clean: results
	@make tests/mhd_clean -C ../..
