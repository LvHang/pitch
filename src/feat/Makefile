

all:

include ../kaldi.mk

TESTFILES = feature-mfcc-test feature-plp-test feature-fbank-test \
         feature-functions-test pitch-functions-test \
         resample-test online-feature-test

OBJFILES = feature-functions.o feature-mfcc.o feature-plp.o feature-fbank.o \
           mel-computations.o wave-reader.o \
           pitch-functions.o resample.o online-feature.o \
		   feature-window.o

LIBNAME = kaldi-feat

ADDLIBS = ../transform/kaldi-transform.a ../util/kaldi-util.a ../thread/kaldi-thread.a \
          ../matrix/kaldi-matrix.a ../base/kaldi-base.a 

include ../makefiles/default_rules.mk

