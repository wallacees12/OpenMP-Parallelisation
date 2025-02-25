CC = gcc
LD = gcc
CFLAGS = -O3 -ftree-vectorize -march=native -ffast-math
LDFLAGS = -lpthread
RM = /bin/rm -f
PROG = galsim
FILE = $(PROG).c

$(PROG): $(FILE)
	$(LD) $(CFLAGS) -o $(PROG) galsim.c $(LDFLAGS) -lm

clean:
	$(RM) $(PROG)