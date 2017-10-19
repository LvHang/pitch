# specific Linux ARM configuration

ifndef DOUBLE_PRECISION
$(error DOUBLE_PRECISION not defined.)
endif

CXXFLAGS = -std=c++11 -I.. $(EXTRA_CXXFLAGS) \
           -Wall -Wno-sign-compare -Wno-unused-local-typedefs \
           -Wno-deprecated-declarations -Winit-self \
           -DKALDI_DOUBLEPRECISION=$(DOUBLE_PRECISION) \
           -DHAVE_EXECINFO_H=1 -DHAVE_CXXABI_H \
           -ftree-vectorize -mfloat-abi=hard -mfpu=neon -pthread \
           -g # -O0 -DKALDI_PARANOID

ifeq ($(PITCH_FLAVOR), dynamic)
CXXFLAGS += -fPIC
endif

# Compiler specific flags
COMPILER = $(shell $(CXX) -v 2>&1)
ifeq ($(findstring clang,$(COMPILER)),clang)
# Suppress annoying clang warnings that are perfectly valid per spec.
CXXFLAGS += -Wno-mismatched-tags
endif

LDFLAGS = $(EXTRA_LDFLAGS) -rdynamic
LDLIBS = $(EXTRA_LDLIBS) -lm -lpthread -ldl
