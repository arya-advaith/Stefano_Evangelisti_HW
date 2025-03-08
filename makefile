# Compiler and flags
CC = gcc
CFLAGS = -Wall -g
LIBS = -llapacke -lm

# Source and object files
OBJECTS = homework.o matrix_solve.o
TARGET = HUCKEL

# Default target
all: $(TARGET)

# Link the program
$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET) $(LIBS)

# Compile source files
homework.o: homework.c matrix_solve.h
	$(CC) $(CFLAGS) -c homework.c

matrix_solve.o: matrix_solve.c matrix_solve.h
	$(CC) $(CFLAGS) -c matrix_solve.c

# Clean up
clean:
	rm -f $(OBJECTS) $(TARGET)

