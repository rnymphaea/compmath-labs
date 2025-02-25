CC=gcc
CFLAGS=-I$(HEADERS) -std=c99 -Wall -Wextra 
MATHFLAG=-lm

BIN=./bin/
HEADERS=../headers

SRC = $(wildcard ./*.c)
OBJ =  $(patsubst ./%.c, $(BIN)%.o, $(SRC))

.PHONY: all clean

all: $(BIN) $(OBJ)
	@ $(CC) $(CFLAGS) $(OBJ) -o main $(MATHFLAG) && ./main

$(BIN):
	@ mkdir $(BIN)

$(BIN)%.o: %.c
	@ $(CC) $(CFLAGS) -c $< -o $@

clean:
	@ rm -rf $(BIN) main
	@ echo "Done!"
