MAINFILE = main
LIBS = ~/workspace/include/toolbox.o ~/workspace/include/sorting.o ~/workspace/include/clock.o
# Intel Fortran Compiler
CFLAGS = -qopenmp -O2 -I/home/$$USER/workspace/include
CC = ~/intel/bin/ifort
# PGI Fortran Compiler
#CFLAGS = -mp -Mbackslash -Minfo -I/home/$$USER/workspace/include
#CC = ~/pgi/linux86-64/2016/bin/pgfortran
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
	$(CC) $(CFLAGS) -c $(MAINFILE).f90 -o Build/$(MAINFILE).o


clean:
	-rm -rf Build
	-rm -f *.mod
