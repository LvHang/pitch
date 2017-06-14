

all:

include ../kaldi.mk

TESTFILES = resample-test

OBJFILES = feature-functions.o mel-computations.o wave-reader.o \
           pitch-functions.o resample.o \
		   feature-window.o

LIBNAME = kaldi-feat

ADDLIBS = ../util/kaldi-util.a \
          ../matrix/kaldi-matrix.a ../base/kaldi-base.a 

include ../makefiles/default_rules.mk

