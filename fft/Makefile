include ../makevars.mk
CFLAGSXX := -openmp $(MKLINC)
LIBS := -openmp $(MKLLIBS)

#########
.PHONY: ../utils
../utils:
	make --directory=$@ all

##########
all: fft_mkl.o fft2D_thrd.o trig.o

fft_mkl.o: fft_mkl.cpp
include fft_mkl.d

fft2D_thrd.o: fft2D_thrd.cpp
include fft2D_thrd.d

trig.o: trig.cpp
include trig.d

testc2c.o: testc2c.cpp
include testc2c.d
testc2c.exe: fft_mkl.o testc2c.o ../utils ../utils/utils.o

testr2c.o: testr2c.cpp
include testr2c.d
testr2c.exe: fft_mkl.o testr2c.o ../utils ../utils/utils.o

test2D.o: test2D.cpp
include test2D.d
test2D.exe: fft_mkl.o test2D.o ../utils 
	$(CPP) -openmp -o $@ fft_mkl.o test2D.o ../utils/utils.o $(MKLTHRD)

test2D_thrd.o: test2D_thrd.cpp
include test2D_thrd.d
test2D_thrd.exe: fft2D_thrd.o fft_mkl.o test2D_thrd.o ../utils ../utils/utils.o


test_trig.o: test_trig.cpp
include test_trig.d
test_trig.exe: test_trig.o trig.o ../utils ../utils/utils.o

time_onemany.o: time_onemany.cpp
include time_onemany.d
time_onemany.exe: fft_mkl.o time_onemany.o ../utils ../utils/utils.o

time_align.o: time_align.cpp
include time_align.d
time_align.exe: fft_mkl.o time_align.o ../utils ../utils/utils.o

time_stride.o: time_stride.cpp
	$(CPP) $(CFLAGS)  -c $<
include time_stride.d
time_stride.exe: fft_mkl.o time_stride.o ../utils ../utils/utils.o

time2D_thrd.o: time2D_thrd.cpp
include time2D_thrd.d
time2D_thrd.exe: fft_mkl.o fft2D_thrd.o time2D_thrd.o ../utils ../utils/utils.o 

time2D_mkl.o: time2D_mkl.cpp
include time2D_mkl.d
time2D_mkl.exe: fft_mkl.o time2D_mkl.o ../utils ../utils/utils.o
	$(CPP) -openmp -o $@ fft_mkl.o time2D_mkl.o ../utils/utils.o $(MKLTHRD)

time_trig.o: time_trig.cpp
include time_trig.d
time_trig.exe: time_trig.o trig.o ../utils ../utils/utils.o
