include ../makevars.mk
CFLAGSXX := $(MKLINC) #-openmp #for test_bvpi
LIBS := $(MKLLIBS) #-openmp #for test_bvpi
MODS := ../utils ../trig ../banded ../pyplot
#########
.PHONY: $(MODS)
$(MODS):
	make --directory=$@ all

##########

bvpi.o: bvpi.cpp
include bvpi.d

test_bvpi.o: test_bvpi.cpp
include test_bvpi.d

test_bvpi.exe: bvpi.o test_bvpi.o 			\
	../utils ../trig ../banded 			\
	../utils/utils.o 				\
	../trig/trig.o					\
	../banded/tridiag.o 				

bvp4.o: bvp4.cpp
include bvp4.d

test_bvp4.o: test_bvp4.cpp
include test_bvp4.d

test_bvp4.exe: test_bvp4.o bvp4.o bvpi.o	\
	../utils ../trig ../banded ../pyplot	\
	../utils/utils.o			\
	../pyplot/pyplot.o			\
	../trig/trig.o				\
	../banded/tridiag.o			
