include ../Makefile.inc

LDFLAGS += -Wl,-rpath,../lib

OBJS = qdft.o
TARGETS = ../qdft

%.o : %.c
	$(CC) $(CFLAGS) -c $<

all : $(TARGETS)
clean :
	$(RM) *.o
	$(RM) $(TARGETS)
dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c) 

$(TARGETS) : $(OBJS)
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ $(OBJS) $(LINKS)
