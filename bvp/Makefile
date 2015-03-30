include ../makevars.mk
CFLAGSXX := $(MKLINC) #-openmp #for test_bvpi
LIBS := $(MKLLIBS) #-openmp #for test_bvpi
MODS := ../utils ../fft ../banded ../pyplot
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
	../utils ../fft ../banded 			\
	../utils/utils.o 				\
	../fft/trig.o					\
	../banded/tridiag.o 				

bvp4fac.o: bvp4fac.cpp
include bvp4fac.d

bvp4si.o: bvp4si.cpp

test_bvp4.o: test_bvp4.cpp
include test_bvp4.d

test_bvp4.exe: test_bvp4.o bvp4si.o bvp4fac.o bvpi.o	\
	../utils ../fft ../banded ../pyplot	\
	../utils/utils.o			\
	../pyplot/pyplot.o			\
	../fft/trig.o				\
	../banded/tridiag.o ../banded/banded.o

bvp4pg.o: bvp4pg.cpp
