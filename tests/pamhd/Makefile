all:
	@$(MAKE) tests/pamhd_executables -C ../..

t: test
test:
	@$(MAKE) tests/pamhd_tests -C ../..

r: results
results:
	@$(MAKE) tests/pamhd_clean_results -C ../..

c: clean
clean: results
	@$(MAKE) tests/pamhd_clean -C ../..
