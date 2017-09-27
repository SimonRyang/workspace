MAINFILE = LifecycleEntrepreneur
CFLAGS = -qopenmp -O2
CC = ~/intel/bin/ifort
TARGET = Build/prog

all: $(TARGET)
	clear
	./$(TARGET)
	-rm -rf Build
	-rm -f *.mod

$(TARGET): Build/$(MAINFILE).o
	$(CC) -qopenmp -O2 Build/$(MAINFILE).o -o $(TARGET)

Build/$(MAINFILE).o: $(MAINFILE).f90
	mkdir -p Build
	$(CC) $(CFLAGS) -c $(MAINFILE).f90 -o Build/$(MAINFILE).o


clean:
	-rm -rf Build
	-rm -f *.mod
