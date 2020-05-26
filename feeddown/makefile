ROOTLIBDIR      = $(shell root-config --libdir)
ROOTLIBS        = $(shell root-config --libs)
ROOTINCL        = $(shell root-config --cflags)
ROOTLDFLAGS     = $(shell root-config --ldflags)
PLATFORM        = $(shell root-config --platform)
ARCH            = $(shell root-config --arch)

ifeq ($(ARCH),macosx64)
	# mac
	CXX = clang++
	CXXOPT = -Wall -pipe -fPIC -pthread -g 
	# SOFLAGS = -dynamiclib -single_module -undefined dynamic_lookup
	SOFLAGS = -dynamiclib -single_module 
	SOSUFFIX = dylib
else
	# linux
	CXX = g++
	CXXOPT = -Wall -fPIC -g -Wno-unknown-pragmas
	SOFLAGS = -shared
	SOSUFFIX = so
endif

libfeeddown.so: feeddown.h feeddown.cpp
	$(CXX) $(SOFLAGS) $(CXXOPT) feeddown.cpp `root-config --cflags --glibs` -I ${PYTHIA8}/include ${PYTHIA8}/lib/libpythia8.$(SOSUFFIX) $(ROOTLIBDIR)/libEGPythia8.so -o libfeeddown.so

make clean:
	rm -f libfeeddown.so