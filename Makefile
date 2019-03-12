
CC = gcc
RM = rm

SRC = src

OFLAGS = -O3 -march=native
WFLAGS = -Wall -Wextra -Wno-unused-variable -Wno-unused-function -Wno-constant-conversion
LDFLAGS = -lz -lbz2 -llzma
CFLAGS = -std=c99 $(OFLAGS) $(WFLAGS) $(LDFLAGS)


all: tinymasker

tinymasker: $(SRC)/tinymasker.c $(SRC)/dozeu.h
	$(CC) $(CFLAGS) -o tinymasker $(SRC)/tinymasker.c $(SRC)/toml.c $(SRC)/re.c

clean:
	$(RM) -f tinymasker

