# SHELL += -x

CXX = g++
CC = gcc         # used for sph2pipe
# CXX = clang++  # Uncomment these lines
# CC = clang     # to build with Clang.

all: check_required_programs atlas
	@echo -e "\n\n"
	@echo "All done OK."

# make sure check_required_programs runs before anything else:
atlas : | check_required_programs

check_required_programs:
	extras/check_dependencies.sh

clean:

distclean:
	rm -rf ATLAS/

atlas: ATLAS/include/cblas.h

ATLAS/include/cblas.h: | atlas3.8.3.tar.gz
	tar xozf atlas3.8.3.tar.gz ATLAS/include

atlas3.8.3.tar.gz:
	wget -T 10 ftp://ftp.vim.org/vol/2/metalab/distributions/tinycorelinux/2.x/tce/src/libatlas/atlas3.8.3.tar.gz || \
    wget -T 10 http://sourceforge.net/projects/math-atlas/files/Stable/3.8.3/atlas3.8.3.tar.gz || \
	wget --no-check-certificate -T 10 -t 3 http://www.danielpovey.com/files/kaldi/atlas3.8.3.tar.gz
