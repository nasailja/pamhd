# to compile and test pick one makefile from makefiles directory
# or write a new one for your environment
ENVIRONMENT_MAKEFILE = makefiles/macosx_macports

#
# The lines below are not intended to be modified by users
#
CXXFLAGS = -std=c++11 -W -Wall -Wextra -pedantic -O3
CPPFLAGS = -I source
include $(ENVIRONMENT_MAKEFILE)

# filled by project makefiles
EXECUTABLES =
TESTS =
RESULTS =
CLEAN =


all: test


include \
  tests/prettyprint/project_makefile \
  tests/volume_range/project_makefile \
  tests/muparserx/project_makefile \
  tests/mhd/project_makefile \
  tests/program_options/project_makefile \
  tests/boundaries/project_makefile


%.tst: %.exe
	@echo "RUN "$< && ./$< && echo PASS && touch $@


t: test
test: $(EXECUTABLES) $(TESTS)
	@echo && echo "All tests passed."

r: results
results: $(RESULTS)

c: clean
clean: results $(CLEAN)
