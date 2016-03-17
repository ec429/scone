CC = gcc
CFLAGS := -Wall -Wextra -Werror --pedantic --std=gnu11 -g
LIBS := -lm
OBJS := vectors.o
INCS := $(OBJS:.o=.h) planets.h

all: scone

scone: scone.c $(OBJS) $(INCS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) $(LIBS) $(OBJS) -o $@

%.o: %.c %.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

planets.h: cfg.py kop.py rssk.cfg
	./kop.py > planets.h
