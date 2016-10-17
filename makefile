CC=gcc
CFLAGS=-g -Wall -c -I/home/nnvv/local/include
LDFLAGS=-lm -L/home/nnvv/local/lib -lgsl -lgslcblas 
SOURCES=allvars.c functions.c main.c 
OBJECTS=$(SOURCES:%.c=%.o)

all: main

main: $(OBJECTS) 
	$(CC)  $(OBJECTS) -o $@.x $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o *.x *~
