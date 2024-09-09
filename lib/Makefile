include ../Makefile.inc

CFLAGS += -fPIC
LDFLAGS += -shared

OBJS = areig.o
TARGETS = libqdft.so

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
	$(CXX) $(LDLIBS) $(LDFLAGS) -o $@ $(OBJS) $(LINKS)
