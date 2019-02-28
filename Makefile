
CC = gcc
RM = rm

SRC = src

OFLAGS = -O3 -march=native
WFLAGS = -Wall -Wextra -Wno-unused-variable -Wno-unused-function -Wno-constant-conversion
CFLAGS = -std=c99 $(OFLAGS) $(WFLAGS)


all: tinymasker

tinymasker: $(SRC)/tinymasker.c $(SRC)/dozeu.h
	$(CC) $(CFLAGS) -o tinymasker $(SRC)/tinymasker.c

clean:
	$(RM) -f tinymasker

