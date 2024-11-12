include ../Makefile.inc

OBJS = fci.o soft.o
TARGETS = ../fci ../soft ../lib/libsoft.so

%.o : %.c
	$(CC) $(CFLAGS) -c $<

all : $(TARGETS)
clean :
	$(RM) *.o
	$(RM) $(TARGETS)
dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c) 

../fci ../soft : ../% : %.o
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ $< $(LINKS)
../lib/libsoft.so : soft.o
	$(CC) $(LDLIBS) $(LDFLAGS) -shared -o $@ $< $(LINKS)
