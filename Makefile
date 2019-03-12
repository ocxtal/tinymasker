
CC = gcc
RM = rm

SRC = src

OFLAGS = -O3
WFLAGS = -Wall -Wextra -Wno-unused-variable -Wno-unused-function -Wno-constant-conversion -Wno-implicit-fallthrough
LDFLAGS = -lz -lbz2 -llzma -lpthread
CFLAGS = -std=c99 -march=native $(OFLAGS) $(WFLAGS) $(LDFLAGS)
DEBUG_CFLAGS = -g -std=c99 -march=native -DDEBUG $(WFLAGS) $(LDFLAGS)


all: tinymasker

tinymasker: $(SRC)/tinymasker.c $(SRC)/dozeu.h
	$(CC) $(CFLAGS) -o tinymasker $(SRC)/tinymasker.c $(SRC)/toml.c $(SRC)/re.c

debug: $(SRC)/tinymasker.c $(SRC)/dozeu.h
	$(CC) $(DEBUG_CFLAGS) -o tinymasker.debug $(SRC)/tinymasker.c $(SRC)/toml.c $(SRC)/re.c

clean:
	$(RM) -f tinymasker tinymasker.debug

