all:
	@$(MAKE) tests/particle_executables -C ../..

t: test
test:
	@$(MAKE) tests/particle_tests -C ../..

r: results
results:
	@$(MAKE) tests/particle_clean_results -C ../..

c: clean
clean: results
	@$(MAKE) tests/particle_clean -C ../..
