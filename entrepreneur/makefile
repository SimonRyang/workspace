
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
	$(CC) $(CFLAGS) $(LIBS) Build/main.o -o $(TARGET)

Build/$(MAINFILE).o: $(MAINFILE).f90
	mkdir -p Build
	$(CC) $(CFLAGS) -I/home/$$USER/repo_projects/include -c globals.f90 -o Build/globals.o
	$(CC) $(CFLAGS) -I/home/$$USER/repo_projects/include -c main.f90 -o Build/main.o

clean:
	-rm -rf Build
	-rm -f *.mod
