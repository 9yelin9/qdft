include ../Makefile.inc

OBJS = areig.o laeig.o util.o
TARGETS = ../lib/libqdft.a

%.o : %.c
	$(CC) $(CFLAGS) -c $<
%.o : %.cpp
	$(CXX) $(CFLAGS) -c $<

all : $(TARGETS)
clean :
	$(RM) *.o
	$(RM) $(TARGETS)
dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c) 

$(TARGETS) : $(OBJS)
	$(AR) $@ $^
