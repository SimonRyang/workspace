MAINFILE = main
LIBS = toolbox.o sorting.o clock.o
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
	$(CC) $(CFLAGS) -o $(TARGET) Build/$(MAINFILE).o $(LIBS)

Build/$(MAINFILE).o: $(MAINFILE).f90
	mkdir -p Build
	$(CC) $(CFLAGS) -c $(MAINFILE).f90 -o Build/$(MAINFILE).o

Build/toolbox.0: toolbox.f90
	$(CC) $(CFLAGS) -c toolbox.f90 -o Build/toolbox.o

Build/toolbox.0: sorting.f90
	$(CC) $(CFLAGS) -c sorting.f90 -o Build/sorting.o

Build/toolbox.0: clock.f90
	$(CC) $(CFLAGS) -c clock.f90 -o Build/clock.o

clean:
	-rm -rf Build
	-rm -f *.mod
