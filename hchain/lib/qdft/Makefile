ESPRESSO_DIR = $(HOME)/espresso
CONDA_DIR = /APP/enhpc/miniconda3

include $(ESPRESSO_DIR)/make.inc

CFLAGS += -I$(CONDA_DIR)/include/python3.12
F90FLAGS += -I$(ESPRESSO_DIR)/UtilXlib -I$(ESPRESSO_DIR)/Modules
#LDFLAGS += -L$(CONDA_DIR)/lib/python3.12/config-3.12-x86_64-linux-gnu -L$(CONDA_DIR)/lib -Wl,-rpath=$(CONDA_DIR)/lib 
#LINKS += -lpython3.12 -lpthread -ldl -lutil -lm 

OBJ = vqe.o qdft.o
TARGET = libqdft.a

.PHONY: all clean

all: $(TARGET)
clean:
	rm -rf $(OBJ) $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) -c $<
%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

libqdft.a: $(OBJ)
	ar rv $@ $^
