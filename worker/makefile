MAINFILE = main
LIBS = ~/repo_projects/include/toolbox.o ~/repo_projects/include/sorting.o ~/repo_projects/include/clock.o
CFLAGS = -O2 -I/home/$$USER/repo_projects/include
CC = ~/intel/bin/ifort
TARGET = Build/prog

all: $(TARGET)
	clear
	./$(TARGET)
	-rm -rf Build
	-rm -f *.mod

$(TARGET): Build/$(MAINFILE).o
	$(CC) -O2 -qopenmp $(LIBS) Build/$(MAINFILE).o -o $(TARGET)

Build/$(MAINFILE).o: $(MAINFILE).f90
	mkdir -p Build
	$(CC) -qopenmp $(CFLAGS) -c $(MAINFILE).f90 -o Build/$(MAINFILE).o


clean:
	-rm -rf Build
	-rm -f *.mod
