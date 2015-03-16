#########
CPP 	:= icpc
CFLAGS 	:= -O3 -prec-div -no-ftz -restrict -Wshadow -fno-inline-functions
MPIINC 	:= `mpiCC -showme:compile`
MKLINC 	:= -I $$MKL_INC


MPILIBS := `mpiCC -showme:link`
MKLLIBS := -L$$MKL_LINK -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
MKLTHRD := -L$$MKL_LINK -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
#MKLINC := -mkl=sequential
#MKLLIBS := -mkl=sequential
#MKLTHRD := -mkl=threaded


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
.PHONY: clean cleanx cleanxx
clean:
	rm *.o; rm *.exe; 

cleanx:
	rm pbs*.*; 

cleanxx:
	rm *.o; rm *.exe; rm pbs*.*; rm *.d; rm DBG/outP*;