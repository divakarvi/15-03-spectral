include ../makevars.mk
CFLAGSXX := -openmp $(MKLINC)
LIBS := -openmp $(MKLLIBS)
MODS := ../utils ../banded ../fft ../globals ../pyplot ../diff
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
	../utils ../fft ../banded ../globals		\
	../utils/utils.o 				\
	../fft/trig.o					\
	../banded/tridiag.o 				\
	../globals/globals.o 				\
	../globals/gl_timer.o

bvp.o: bvp.cpp
include bvp.d

test_bvp.exe: bvp.o bvpi.o test_bvp.o	 			\
	../utils ../fft ../banded ../globals ../pyplot ../diff	\
	../utils/utils.o 					\
	../fft/trig.o						\
	../banded/tridiag.o 					\
	../banded/banded.o 					\
	../globals/globals.o 					\
	../globals/gl_timer.o					\
	../pyplot/pyplot.o					\
	../diff/diff_fd.o					\
	../diff/diff_trig.o					\
	../diff/fdweights.o					

error_graph.exe: bvp.o bvpi.o error_graph.o	 		\
	../utils ../fft ../banded ../globals ../pyplot ../diff	\
	../utils/utils.o 					\
	../fft/trig.o						\
	../banded/tridiag.o 					\
	../banded/banded.o 					\
	../globals/globals.o 					\
	../globals/gl_timer.o					\
	../pyplot/pyplot.o					\
	../diff/diff_trig.o					\
	../diff/diff_fd.o					\
	../diff/fdweights.o					
