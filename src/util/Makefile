

all:


include ../pitch.mk

TESTFILES = text-utils-test kaldi-io-test parse-options-test \
    kaldi-table-test

OBJFILES = text-utils.o kaldi-io.o \
        kaldi-holder.o  kaldi-table.o parse-options.o simple-io-funcs.o

LIBNAME = pitch-util

ADDLIBS = ../matrix/pitch-matrix.a \
          ../base/pitch-base.a 

include ../makefiles/default_rules.mk
