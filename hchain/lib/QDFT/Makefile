CC = gcc

CONDA_DIR = /APP/enhpc/miniconda3

CFLAGS = -gdwarf-2 -O2 -Wall -mcmodel=medium -I$(CONDA_DIR)/include/python3.12
LDFLAGS = -L$(CONDA_DIR)/lib/python3.12/config-3.12-x86_64-linux-gnu -L$(CONDA_DIR)/lib -Wl,-rpath=$(CONDA_DIR)/lib
LINKS = -lpython3.12 -lpthread -ldl -lutil -lm

OBJ = qdft_wrapper.o
TARGET = libqdft.a #test

.PHONY: all clean

all: $(TARGET)
clean:
	rm -rf $(OBJ) $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

libqdft.a: $(OBJ)
	ar rv $@ $<
#test: $(OBJ)
#	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ $< $(LINKS)
