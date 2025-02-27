CC = gcc
LD = gcc
CFLAGS = -O3 -ftree-vectorize -march=native -ffast-math
LDFLAGS = -lpthread -lm
RM = /bin/rm -f
PROG_BALANCED = galsim_balanced
PROG_UNBALANCED = galsim_unbalanced
FILE_BALANCED = galsim_balanced.c
FILE_UNBALANCED = galsim_unbalanced.c

# Default target (build both balanced and unbalanced versions)
all: $(PROG_BALANCED) $(PROG_UNBALANCED)

# Rule for building galsim_balanced
$(PROG_BALANCED): $(FILE_BALANCED)
	$(LD) $(CFLAGS) -o $(PROG_BALANCED) $(FILE_BALANCED) $(LDFLAGS)

# Rule for building galsim_unbalanced
$(PROG_UNBALANCED): $(FILE_UNBALANCED)
	$(LD) $(CFLAGS) -o $(PROG_UNBALANCED) $(FILE_UNBALANCED) $(LDFLAGS)

# Clean up build artifacts
clean:
	$(RM) $(PROG_BALANCED) $(PROG_UNBALANCED)