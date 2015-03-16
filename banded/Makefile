include ../makevars.mk
CFLAGSXX := $(MKLINC)
LIBS := $(MKLLIBS)

#########
.PHONY: ../utils
../utils:
	make --directory=$@ all

##########
.PHONY: all

all: tridiag.o banded.o

tridiag.o: tridiag.cpp
include tridiag.d

banded.o: banded.cpp
include banded.d

test_tridiag.o: test_tridiag.cpp
include test_tridiag.d

test_tridiag.exe: tridiag.o test_tridiag.o ../utils ../utils/utils.o

test_banded.o: test_banded.cpp
include test_banded.d

test_banded.exe: banded.o test_banded.o ../utils ../utils/utils.o

time_tri.o: time_tri.cpp
include time_tri.d

time_tri.exe: tridiag.o time_tri.o ../utils ../utils/utils.o

time_banded.o: time_banded.cpp
include time_banded.d

time_banded.exe: banded.o time_banded.o ../utils ../utils/utils.o
