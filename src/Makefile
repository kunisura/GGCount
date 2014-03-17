# Computing the Number of Paths in a Grid Graph
# Hiroaki Iwashita <iwashita@erato.ist.hokudai.ac.jp>
# Copyright (c) 2013 ERATO MINATO Project
# $Id$

TARGET        = ggcount
SINGLE_TARGET = $(TARGET:%=%_single)
MULTI_TARGET  = $(TARGET:%=%_multi)
DEBUG_TARGET  = $(TARGET:%=%_debug)

CPPFLAGS      = -std=c++11
CXXFLAGS      = -Wall -fmessage-length=0

.PONY:  all single multi debug clean

all:    single multi
single: $(SINGLE_TARGET)
multi:  $(MULTI_TARGET)
debug:  $(DEBUG_TARGET)

clean:
	$(RM) $(SINGLE_TARGET) $(MULTI_TARGET) $(DEBUG_TARGET) *.o

%_single: %.cpp
	$(LINK.cc) -DNDEBUG -O3 $^ $(LOADLIBES) $(LDLIBS) -o $@

%_multi: %.cpp
	$(LINK.cc) -DNDEBUG -O3 -fopenmp $^ $(LOADLIBES) $(LDLIBS) -o $@

%_debug: %.cpp
	$(LINK.cc) -g $^ $(LOADLIBES) $(LDLIBS) -o $@

%: %.o
	$(CXX) $(LDFLAGS) $(TARGET_ARCH) $^ $(LDLIBS) -o $@
