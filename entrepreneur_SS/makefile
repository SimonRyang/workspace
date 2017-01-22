MAINFILE = main
LIBS = globals.o ~/workspace/include/toolbox.o ~/workspace/include/sorting.o ~/workspace/include/clock.o

# Intel Fortran Compiler
CFLAGS = -qopenmp -O2
CC = ~/intel/bin/ifort

# PGI Fortran Compiler
#CFLAGS = -mp -O2 -Mbackslash
#CC = ~/pgi/linux86-64/2016/bin/pgfortran

# GFortran Compiler
#CFLAGS = -fopenmp -O2
#CC = gfortran

TARGET = prog

all: $(TARGET)
	clear
	./$(TARGET)
	-rm -f *.mod

$(TARGET): $(MAINFILE).o $(LIBS)
	$(CC) $(CFLAGS) $(MAINFILE).o $(LIBS) -o $(TARGET) 

globals.mod: ~/workspace/include/toolbox.mod ~/workspace/include/sorting.mod globals.f90
	$(CC) $(CFLAGS) -c globals.f90

~/workspace/include/toolbox.mod: ~/workspace/include/toolbox.f90
	$(CC) $(CFLAGS) -c ~/workspace/include/toolbox.f90

~/workspace/include/sorting.mod: ~/workspace/include/toolbox.mod ~/workspace/include/sorting.f90
	$(CC) $(CFLAGS) -c ~/workspace/include/sorting.f90

~/workspace/include/clock.mod: ~/workspace/include/clock.f90
	$(CC) $(CFLAGS) -c ~/workspace/include/clock.f90

$(MAINFILE).o: globals.mod ~/workspace/include/toolbox.mod ~/workspace/include/sorting.mod ~/workspace/include/clock.mod $(MAINFILE).f90
	$(CC) $(CFLAGS) -c $(MAINFILE).f90

#clean:
#	-rm -f *.mod
#	-rm -f *.o
