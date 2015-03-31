include ../makevars.mk
CFLAGSXX := -openmp $(MKLINC)
LIBS := -openmp $(MKLLIBS) -lifcore
MODS := ../utils ../pyplot ../fft
########
.PHONY: $(MODS)
$(MODS):
	make --directory=$@ all

########
all: legcheb.o legendre.o 

legcheb.o: legcheb.f
	ifort -O3 -c $<

legendre.o: legendre.cpp

test_legendre.o: test_legendre.cpp

test_legendre.exe: test_legendre.o legendre.o legcheb.o		\
	../utils ../utils/utils.o				\
	../pyplot ../pyplot/pyplot.o				\
	../fft ../fft/trig.o

time_legendre.o: time_legendre.cpp

time_legendre.exe: time_legendre.o legendre.o legcheb.o		\
	../utils ../utils/utils.o				\
	../pyplot ../pyplot/pyplot.o				\
	../fft ../fft/trig.o
