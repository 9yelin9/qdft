include ../Makefile.inc

OBJS = qdft.o test.o
TARGETS = ../qdft ../test

%.o : %.c
	$(CC) $(CFLAGS) -c $<

all : $(TARGETS)
clean :
	$(RM) *.o
	$(RM) $(TARGETS)
dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c) 

$(TARGETS) : ../% : %.o
	$(CXX) $(LDLIBS) $(LDFLAGS) -o $@ $< $(LINKS)
