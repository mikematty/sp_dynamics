############################################################
#               Stochastic LLG project                     #
#                    O. Hovorka                            #
#                     10/8/2011                            #
############################################################
# home directory
BASE = /Users/ohovorka/MyWork/Programming/Cpp/Random_numbers

############################################################

# object file 
OBJMAIN = LLGlangevin3D.o \
          mersenne.o \
          userintf.o \
          stoc1.o

############################################################


# executable file
EXE = run

# macro for name of my compiler and linker
CC = g++ 
LD = $(CC) 

# macro for compilation options:
# directory for include files
CFLAG = -O3 -msse3 -funroll-loops -funroll-all-loops \
            -I Random_numbers/randomc \
            -I Random_numbers/stocc

# include libraries (.so,.a,etc) here...
LDFLAG = -lfftw3 -lm

whatever: $(OBJMAIN)
	$(LD) $(LDFLAG) -o $(EXE) $(OBJMAIN)
LLGlangevin3D.o:
	$(CC) $(CFLAG) -c $*.cpp -o $*.o
userintf.o:
	$(CC) $(CFLAG) -c Random_numbers/randomc/$*.cpp -o $*.o
mersenne.o:
	$(CC) $(CFLAG) -c Random_numbers/randomc/$*.cpp -o $*.o
stoc1.o:
	$(CC) $(CFLAG) -c Random_numbers/stocc/$*.cpp -o $*.o
clean:	
	@rm -f $(OBJMAIN) $(EXE)
