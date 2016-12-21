MAINFILE = main
LIBS = globals.mod toolbox.mod sorting.mod clock.mod
# Intel Fortran Compiler
CFLAGS = -qopenmp -O2
CC = ~/intel/bin/ifort
# PGI Fortran Compiler
#CFLAGS = -mp -Mbackslash
#CC = ~/pgi/linux86-64/2016/bin/pgfortran
TARGET = prog

all: $(TARGET)
	clear
	./$(TARGET)
	-rm -f *.mod

$(TARGET): $(MAINFILE).o $(LIBS)
	$(CC) $(CFLAGS) $(MAINFILE).o $(LIBS) -o $(TARGET) 

globals.mod: sorting.mod globals.f90
	$(CC) $(CFLAGS) -c globals.f90 -o globals.mod

toolbox.mod: toolbox.f90
	$(CC) $(CFLAGS) -c toolbox.f90 -o toolbox.mod

sorting.mod: sorting.f90
	$(CC) $(CFLAGS) -c sorting.f90 -o sorting.mod

clock.mod: clock.f90
	$(CC) $(CFLAGS) -c clock.f90 -o clock.mod

$(MAINFILE).o: globals.mod toolbox.mod sorting.mod clock.mod $(MAINFILE).f90
	$(CC) $(CFLAGS) -c $(MAINFILE).f90

clean:
	-rm -rf Build
	-rm -f *.mod
