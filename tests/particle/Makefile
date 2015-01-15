all:
	@make tests/particle_executables -C ../..

t: test
test:
	@make tests/particle_tests -C ../..

r: results
results:
	@make tests/particle_clean_results -C ../..

c: clean
clean: results
	@make tests/particle_clean -C ../..
