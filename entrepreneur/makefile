MAINFILE = main
LIBS = toolbox.o sorting.o clock.o

# Intel Fortran Compiler
CFLAGS = -qopenmp -O2 #-check all
CC = ~/intel/bin/ifort

TARGET = prog

all: $(TARGET)
	clear
	./$(TARGET)

$(TARGET): $(MAINFILE).o $(LIBS)
	$(CC) $(CFLAGS) $(LIBS) $(MAINFILE).o  -o $(TARGET)

$(MAINFILE).o: $(MAINFILE).f90
	$(CC) $(CFLAGS) -c $(MAINFILE).f90

clean:
	-rm -f main.o
