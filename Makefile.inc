CC = gcc
CXX = g++
AR = ar rv
RM = rm -rf

USR_DIR = /home/yerin/.local
ARPACK_DIR = $(USR_DIR)/arpackpp-2.4.0

CFLAGS = -g -O2 -Wall -mcmodel=medium -I../include -I$(USR_DIR)/include -I$(ARPACK_DIR)/include -fopenmp
LDFLAGS = -L../lib -L$(USR_DIR)/lib -L$(USR_DIR)/lib64 -fopenmp 
LINKS = -lqdft -lz -lm -lopenblas -llapack -larpack -lsuperlu -lgfortran

.PHONY: all clean dep
