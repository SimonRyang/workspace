MAINFILE = main
LIBS = globals.o toolbox.o sorting.o clock.o

# Intel Fortran Compiler
CFLAGS = -qopenmp -O2 #-check all
CC = ~/intel/bin/ifort

TARGET = prog

all: $(TARGET)
	clear
	./$(TARGET)
	-rm -f globals.mod

$(TARGET): $(MAINFILE).o $(LIBS)
	$(CC) $(CFLAGS) $(MAINFILE).o $(LIBS) -o $(TARGET)

globals.mod: toolbox.mod sorting.mod globals.f90
	$(CC) $(CFLAGS) -c globals.f90

toolbox.mod: toolbox.f90
	$(CC) $(CFLAGS) -c toolbox.f90

sorting.mod: toolbox.mod sorting.f90
	$(CC) $(CFLAGS) -c sorting.f90

clock.mod: clock.f90
	$(CC) $(CFLAGS) -c clock.f90

$(MAINFILE).o: globals.mod toolbox.mod sorting.mod clock.mod $(MAINFILE).f90
	$(CC) $(CFLAGS) -c $(MAINFILE).f90

clean:
	-rm -f globals.mod
	-rm -f globasl.o main.o
