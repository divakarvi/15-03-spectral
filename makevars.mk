#########
CPP 	:= icpc
CFLAGS 	:= -O3 -prec-div -no-ftz -restrict -Wshadow -fno-inline-functions

#MKLINC 	:= -I $$MKL_INC
#MKLLIBS := -L$$MKL_LINK -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
#MKLTHRD := -L$$MKL_LINK -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm

MKLINC := -mkl=sequential
MKLLIBS := -mkl=sequential
MKLTHRD := -mkl=threaded


#########
.SUFFIXES: .cpp .c .o .exe .s .d

%.o: %.cpp
	$(CPP) $(CFLAGS) $(CFLAGSXX) -c $<
%.d: %.cpp
	@set -e; rm -f $@; \
	$(CPP) -M $(CFLAGS) $(CFLAGSXX) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
%.s: %.cpp 
	$(CPP) $(CFLAGS) $(CFLAGSXX) -S $< 
%.exe: %.o 
	$(CPP) -o $@ $(filter %.o,$^) $(LIBS)

########
.PHONY: clean cleanxx
clean:
	rm *.o; rm *.exe; 
cleanxx:
	rm *.o; rm *.exe; rm pbs*.*; rm *.d; rm DBG/outP*;