

all:

include ../pitch.mk

TESTFILES = resample-test

OBJFILES = feature-functions.o wave-reader.o \
           pitch-functions.o resample.o \
		   feature-window.o

LIBNAME = pitch-feat

ADDLIBS = ../util/pitch-util.a \
          ../matrix/pitch-matrix.a ../base/pitch-base.a 

include ../makefiles/default_rules.mk

