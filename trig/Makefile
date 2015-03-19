include ../makevars.mk
CFLAGSXX := -openmp $(MKLINC)
LIBS := -openmp $(MKLLIBS)

#########
.PHONY: ../utils
../utils:
	make --directory=$@ all

##########
all: trig.o

trig.o: trig.cpp
include trig.d

test_trig.o: test_trig.cpp
include test_trig.d
test_trig.exe: test_trig.o trig.o ../utils ../utils/utils.o

time_trig.o: time_trig.cpp
include time_trig.d
time_trig.exe: time_trig.o trig.o ../utils ../utils/utils.o
