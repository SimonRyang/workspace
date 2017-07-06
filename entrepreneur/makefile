MAINFILE = main
LIBS = /home/flh71wr/workspace/include/toolbox.o /home/flh71wr/workspace/include/sorting.o /home/flh71wr/workspace/include/clock.o
CFLAGS = -qopenmp -O2 -I/home/flh71wr/workspace/include
CC = ~/intel/bin/ifort
TARGET = Build/prog

all: $(TARGET)
	clear
	./$(TARGET)
	-rm -rf Build
	-rm -f *.mod

$(TARGET): Build/$(MAINFILE).o
	$(CC) -qopenmp -O2 $(LIBS) Build/$(MAINFILE).o -o $(TARGET)

Build/$(MAINFILE).o: $(MAINFILE).f90
	mkdir -p Build
	$(CC) $(CFLAGS) -c $(MAINFILE).f90 -o Build/$(MAINFILE).o


clean:
	-rm -rf Build
	-rm -f *.mod
