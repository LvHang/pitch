# Cygwin configuration

ifndef DOUBLE_PRECISION
$(error DOUBLE_PRECISION not defined.)
endif

CXXFLAGS = -std=c++11 -I.. $(EXTRA_CXXFLAGS) \
           -Wall -Wno-sign-compare -Wno-unused-local-typedefs \
           -Wno-deprecated-declarations -Winit-self \
           -DKALDI_DOUBLEPRECISION=$(DOUBLE_PRECISION) \
           -msse -msse2 \
           -g # -O0 -DKALDI_PARANOID

ifeq ($(KALDI_FLAVOR), dynamic)
CXXFLAGS += -fPIC
endif

LDFLAGS = $(EXTRA_LDFLAGS) -g --enable-auto-import
LDLIBS = $(EXTRA_LDLIBS) -lm -lpthread -ldl
