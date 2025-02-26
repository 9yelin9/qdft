CC = gcc
CXX = g++
AR = ar rv
RM = rm -rf

LOCAL_DIR = $(HOME)/.local
ARPACK_DIR = $(HOME)/.arpackpp-2.4.0

CFLAGS = -g -O2 -Wall -mcmodel=medium -fPIC -I../include -I$(LOCAL_DIR)/include -I$(ARPACK_DIR)/include -fopenmp
LDFLAGS = -L../lib -L$(LOCAL_DIR)/lib -L$(LOCAL_DIR)/lib64 -fopenmp 
LINKS = -lhub1d -lz -lm -lopenblas -llapack -larpack -lsuperlu -lgfortran -lgsl

.PHONY: all clean dep
