include ../Makefile.inc

OBJS = soft.o fci.o test.o
TARGETS = ../soft ../fci ../test

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
