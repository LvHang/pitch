# This is the top-level Makefile for standalone pitch extractor.
# Also see pitch.mk which supplies options and some rules
# used by the Makefiles in the subdirectories.

SHELL := /bin/bash


SUBDIRS = base matrix util feat featbin

MEMTESTDIRS = base matrix util feat featbin

SUBDIRS_LIB = $(filter-out %bin, $(SUBDIRS))

PITCH_SONAME ?= libpitch.so

include pitch.mk

# Reset the default goal, so that the all target will become default
.DEFAULT_GOAL :=
all: 
	$(MAKE) checkversion 
	$(MAKE) pitch.mk 
	$(MAKE) mklibdir
	$(MAKE) subdirs
	-echo Done

# In ./configure file, it set a variable to denote the version
# The "CONFIGURE_VERSION" variable is set in pitch.mk
.PHONY: checkversion
checkversion:
ifeq ($(shell ./configure --version),$(CONFIGURE_VERSION))
	@echo "The version of configure script matches pitch.mk version. Good."
else
	@echo ""
	@echo "The pitch.mk file was generated using a different version of configure script. Please rerun the configure again"
	@test -f ./pitch.mk && echo  "Hint: Previous configure command line: " && head -n 2 ./pitch.mk | grep configure | sed 's/^# *//g'
	@echo ""
	@false
endif


pitch.mk:
	@[ -f pitch.mk ] || { echo "pitch.mk does not exist; you have to run ./configure"; exit 1; }

# The "PITCHLIBDIR" variable is set in pitch.mk
# If the directory is not exist, make a dir
mklibdir:
	test -d $(PITCHLIBDIR) || mkdir $(PITCHLIBDIR)

subdirs: $(SUBDIRS)

#I don't want to call rm -rf
#If pitch lib dir is exist, remove the *.so files firstly, and then remove the dir.
rmlibdir:
ifneq ($(PITCHLIBDIR), )
	-rm $(PITCHLIBDIR)/*{.so,.a,.o}
	-rmdir $(PITCHLIBDIR)
else
	@true
endif

#SUBDIRS depends on "mklibdir", so we need to make "mklibdir" first.
#$@ means target file. So "make -C" go into the subdir.
.PHONY: $(SUBDIRS)
$(SUBDIRS) : mklibdir
	$(MAKE) -C $@

#PITCH_FLAVOR is set in pitch.mk. We will build the biglib depend on different system.
#CXX is set in pitch.mk and equals to g++ as default.
#The rule make the *.dylib files in subdirs are merged to libkaldi.dylib
#for linux system, the */kaldi-*.a, libm.so, libpthread.so, libdl.so will be combine
biglib: $(SUBDIRS_LIB)
ifeq ($(PITCH_FLAVOR), dynamic)
ifeq ($(shell uname), Darwin)
	$(CXX) -dynamiclib -o $(PITCHLIBDIR)/libkaldi.dylib -install_name @rpath/libkaldi.dylib -framework Accelerate $(LDFLAGS) $(SUBDIRS_LIB:=/*.dylib)
else
ifeq ($(shell uname), Linux)
	#$(warning the following command will probably fail, in that case add -fPIC to your CXXFLAGS and remake all)
	$(CXX) -shared -o $(PITCHLIBDIR)/$(PITCH_SONAME) -Wl,-soname=$(PITCH_SONAME),--whole-archive  $(SUBDIRS_LIB:=/pitch-*.a) $(LDLIBS) -Wl,--no-whole-archive
else
	$(error Dynamic libraries not supported on this platform. Run configure with --static flag. )
endif
endif
endif

#Firstly, do rmlibdir
#Then, go into the subdir and do make clean.
clean: rmlibdir
	-for x in $(SUBDIRS); do $(MAKE) -C $$x clean; done

#distclean means rm the .depend.mk
distclean: clean
	-for x in $(SUBDIRS); do $(MAKE) -C $$x distclean; done

#prerequisite is subdir/test
test: $(addsuffix /test, $(SUBDIRS))

# Define an implicit rule, expands to e.g.:
#  base/test: base
#     $(MAKE) -C base test
#prerequisite is the dir and mklibdir(make library dir)
%/test: % mklibdir
	$(MAKE) -C $< test

valgrind:
	-for x in $(MEMTESTDIRS); do $(MAKE) -C $$x valgrind || { echo "valgrind on $$x failed"; exit 1; }; done

#It is similar with test rules.
depend: $(addsuffix /depend, $(SUBDIRS))

base/.depend.mk:
	$(MAKE) depend

%/depend:
	$(MAKE) -C $(dir $@) depend

### Dependency list ###
# this is necessary for correct parallel compilation
#1)The tools depend on all the libraries
featbin: base matrix util feat

#2)The libraries have inter-dependencies
base: base/.depend.mk
matrix: base
util: base matrix
feat: base matrix util

