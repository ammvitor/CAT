#!/bin/bash
CC=g++
CFLAGS=-c -w -std=c++11 -m64 -pipe -O2 -Wall -W
SOURCES=src/cosolv_pdb.cpp src/receptor_pdb.cpp src/main.cpp
LDFLAGS        = -m64 -Wl,-O1
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE= CAT

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm $(OBJECTS)

