CC = gcc

CONDA_DIR = /APP/enhpc/miniconda3

CFLAGS = -g -O2 -Wall -mcmodel=medium -fPIC -I$(CONDA_DIR)/include/python3.12
LDFLAGS = -L$(CONDA_DIR)/lib/python3.12/config-3.12-x86_64-linux-gnu -L$(CONDA_DIR)/lib -Wl,-rpath=$(CONDA_DIR)/lib
LINKS = -lpython3.12 -lpthread -ldl -lutil -lm

OBJ = qdft_wrapper.o
TARGET = libqdft.so

.PHONY: all clean

all: $(TARGET)
clean:
	rm -rf *.o $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -shared -o $@ $< $(LDFLAGS) $(LINKS)
