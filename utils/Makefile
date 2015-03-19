include ../makevars.mk

#########
all: utils.o Table.o
utils.o: utils.cpp
include utils.d
Table.o: Table.cpp
include Table.d
test_stat.exe: test_stat.o
test_stat.o: test_stat.cpp
test_utils.exe: test_utils.o utils.o
test_utils.o: test_utils.cpp
include test_utils.d
include test_stat.d
#########

