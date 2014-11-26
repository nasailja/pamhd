all:
	@$(MAKE) tests/mhd_executables -C ../..

t: test
test:
	@$(MAKE) tests/mhd_tests -C ../..

.PHONY: results
r: results
results:
	@$(MAKE) tests/mhd_clean_results -C ../..

c: clean
clean: results
	@$(MAKE) tests/mhd_clean -C ../..
