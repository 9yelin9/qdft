DIRS = lib main

.PHONY: all clean dep
.SUFFIXES : .c .o

define run_make
	@for d in $(DIRS); \
	do \
		make -C $$d $(1); \
	done
endef

all :
	$(call run_make, )
clean :
	$(call run_make, clean)
dep :
	$(call run_make, dep)
