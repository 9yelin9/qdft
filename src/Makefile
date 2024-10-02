include ../Makefile.inc

OBJS = soft.o fci.o
TARGETS = ../soft ../fci

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
