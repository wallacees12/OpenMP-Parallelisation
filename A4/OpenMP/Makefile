CC = gcc
CFLAGS = -O3  -fopenmp -lm -ffast-math -march=native -ftree-vectorize

TARGET = galsim
SRCS = galsim.c

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)
