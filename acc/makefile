MAINFILE = main
CFLAGS = -mp
CC = /opt/pgi/linux86-64/2016/bin/pgfortran
TARGET = Build/prog

all: $(TARGET)
	clear
	./$(TARGET)
	-rm -rf Build
	-rm -f *.mod

$(TARGET): Build/$(MAINFILE).o
	$(CC) $(CFLAGS) Build/$(MAINFILE).o -o $(TARGET)

Build/$(MAINFILE).o: $(MAINFILE).f90
	mkdir -p Build
	$(CC) $(CFLAGS) -c $(MAINFILE).f90 -o Build/$(MAINFILE).o

clean:
	-rm -rf Build
	-rm -f *.mod
