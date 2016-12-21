MAINFILE = main
LIBS = globals.mod toolbox.o sorting.o clock.o
# Intel Fortran Compiler
CFLAGS = -qopenmp -O2
CC = ~/intel/bin/ifort
# PGI Fortran Compiler
#CFLAGS = -mp -Mbackslash
#CC = ~/pgi/linux86-64/2016/bin/pgfortran
TARGET = Build/prog

all: $(TARGET)
	clear
	./$(TARGET)
	-rm -rf Build
	-rm -f *.mod

$(TARGET): Build/$(MAINFILE).o Build/$(LIBS)
	$(CC) $(CFLAGS) Build/$(MAINFILE).o $(LIBS) -o $(TARGET) 

Build/$(MAINFILE).o: $(MAINFILE).f90
	mkdir -p Build
	$(CC) $(CFLAGS) -c $(MAINFILE).f90 -o Build/$(MAINFILE).o

Build/toolbox.o: toolbox.f90
	$(CC) $(CFLAGS) -c toolbox.f90 -o Build/toolbox.o

Build/sorting.o: sorting.f90
	$(CC) $(CFLAGS) -c sorting.f90 -o Build/sorting.o

Build/clock.o: clock.f90
	$(CC) $(CFLAGS) -c clock.f90 -o Build/clock.o

Build/globals.mod : globals.f90
	$(CC) $(CFLAGS) -c globals.f90 -o Build/globals.mod

clean:
	-rm -rf Build
	-rm -f *.mod
