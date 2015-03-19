include ../makevars.mk
CFLAGSXX := $(MKLINC) #-openmp #for test_bvpi
LIBS := $(MKLLIBS) #-openmp #for test_bvpi
MODS := ../utils ../trig ../banded
#########
.PHONY: $(MODS)
$(MODS):
	make --directory=$@ all

##########
all: bvpi.o bvp.o

bvpi.o: bvpi.cpp
include bvpi.d

test_bvpi.o: test_bvpi.cpp
include test_bvpi.d

test_bvpi.exe: bvpi.o test_bvpi.o 			\
	../utils ../trig ../banded 			\
	../utils/utils.o 				\
	../trig/trig.o					\
	../banded/tridiag.o 				

