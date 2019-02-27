
CC = gcc
RM = rm

SRC = src

OFLAGS = -O3 -march=native
WFLAGS = -Wall -Wextra -Wno-unused-variable -Wno-unused-function -Wno-constant-conversion
CFLAGS = -std=c99 $(OFLAGS) $(WFLAGS)


all: masker

masker: $(SRC)/masker.c $(SRC)/dozeu.h
	$(CC) $(CFLAGS) -o masker $(SRC)/masker.c

clean:
	$(RM) -f masker

