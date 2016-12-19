MAINFILE = main
LIBS = ~/workspace/include/toolbox.f90 ~/workspace/include/sorting.f90 ~/workspace/include/clock.f90
CFLAGS = -O2 -qopenmp -traceback 
CC = ~/intel/bin/ifort
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
	$(CC) $(CFLAGS) $(LIBS) -c $(MAINFILE).f90 -o Build/$(MAINFILE).o

clean:
	-rm -rf Build
	-rm -f *.mod
