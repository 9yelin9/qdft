CC = gcc
CXX = g++
RM = rm -rf

USR_DIR = /home/yerin/.local
ARPACK_DIR = $(USR_DIR)/arpackpp-2.4.0

CFLAGS = -g -O2 -Wall -mcmodel=medium -I../include -I$(USR_DIR)/include -I$(ARPACK_DIR)/include -fopenmp
LDFLAGS = -L$(USR_DIR)/lib -fopenmp 
LINKS = -lz -lm -lopenblas -llapack

.PHONY: all clean dep
