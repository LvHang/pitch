
all:

EXTRA_CXXFLAGS = -Wno-sign-compare

include ../pitch.mk

BINFILES = compute-and-process-kaldi-pitch-feats interpolate-pitch \
		   process-kaldi-pitch-feats compute-kaldi-pitch-feats \
		   process-pitch-feats

OBJFILES =

TESTFILES =

ADDLIBS = ../feat/pitch-feat.a \
		  ../util/pitch-util.a \
          ../matrix/pitch-matrix.a ../base/pitch-base.a 

include ../makefiles/default_rules.mk

