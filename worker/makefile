MAINFILE = main
LIBS = ~/workspace/include/toolbox.f90 ~/workspace/include/sorting.f90 ~/workspace/include/clock.f90
# Intel Fortran Compiler
#CFLAGS = -qopenmp -O2
#CC = ~/intel/bin/ifort
# PGI Fortran Compiler
CFLAGS = -mp -Mbackslash -I/home/$$USER/workspace/include
CC = ~/pgi/linux86-64/2016/bin/pgfortran
TARGET = Build/prog

all: $(TARGET)
	clear
	./$(TARGET)
	-rm -rf Build
	-rm -f *.mod

$(TARGET): Build/$(MAINFILE).o
	$(CC) $(CFLAGS) $(LIBS) Build/$(MAINFILE).o -o $(TARGET)

Build/$(MAINFILE).o: $(MAINFILE).f90
	mkdir -p Build
	$(CC) $(CFLAGS) -I/home/$$USER/workspace/include -c $(MAINFILE).f90 -o Build/$(MAINFILE).o


clean:
	-rm -rf Build
	-rm -f *.mod
