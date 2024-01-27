
DEBUG = 0

#######################################################################################

OPT = -O2

FLAGS = -Wall -Wno-deprecated -I..
#-Wall 

ifeq ($(DEBUG),0)
	FLAGS += $(OPT) -g
	FFLAGS = $(OPT) -g
else
	FLAGS += -g -DDEBUG
	FFLAGS = -g -Wall
endif

CC = g++
FC = gfortran

LD = $(FC)
PROBSRC := $(CURDIR)/PROBLEMS
OBJS := MainGlob.o prog_riB-2.o lsrch_box.o qsortd.o problem.o
LIBS = 

o = .o

.SUFFIXES: .f .f90

all:    ddfsa

ddfsa:	$(OBJS) 
	$(LD) $(FLAGS) -o ddfsa $(OBJS) $(LIBS) -pthread

.f.o:
	$(FC) $(FFLAGS) -c $*.f -o $*.o

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90 -o $*.o

clean:
	rm -f risultato.tex
	rm -f fort.1
	rm -f *.mod
	rm -f *.o
	rm -f ddfsa

