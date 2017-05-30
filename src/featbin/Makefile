
all:
EXTRA_CXXFLAGS = -Wno-sign-compare
include ../kaldi.mk

BINFILES = compute-and-process-kaldi-pitch-feats interpolate-pitch \
		   process-kaldi-pitch-feats compute-kaldi-pitch-feats \
		   process-pitch-feats

OBJFILES =

TESTFILES =

ADDLIBS = ../feat/kaldi-feat.a ../transform/kaldi-transform.a \
		  ../util/kaldi-util.a ../thread/kaldi-thread.a \
          ../matrix/kaldi-matrix.a ../base/kaldi-base.a 

include ../makefiles/default_rules.mk

