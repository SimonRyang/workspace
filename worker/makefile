MAINFILE = main
LIBS = globals.o toolbox.o sorting.o clock.o
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
	$(CC) $(CFLAGS) Build/$(MAINFILE).o $(LIBS) -o $(TARGET) 

$(MAINFILE).o: $(MAINFILE).f90
	mkdir -p Build
	$(CC) $(CFLAGS) -c $(MAINFILE).f90 -o $(MAINFILE).o

toolbox.o: toolbox.f90
	$(CC) $(CFLAGS) -c toolbox.f90 -o toolbox.o

sorting.o: sorting.f90
	$(CC) $(CFLAGS) -c sorting.f90 -o sorting.o

clock.o: clock.f90
	$(CC) $(CFLAGS) -c clock.f90 -o clock.o

globals.o : globals.f90
	$(CC) $(CFLAGS) -c globals.f90 -o globals.o

clean:
	-rm -rf Build
	-rm -f *.mod
