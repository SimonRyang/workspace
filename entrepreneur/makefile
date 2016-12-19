MAINFILE = main
LIBS = ~/repo_projects/include/toolbox.f90 ~/repo_projects/include/sorting.f90 ~/repo_projects/include/clock.f90
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
	$(CC) $(CFLAGS) -I/home/$$USER/repo_projects/include -c $(MAINFILE).f90 -o Build/$(MAINFILE).o

clean:
	-rm -rf Build
	-rm -f *.mod
