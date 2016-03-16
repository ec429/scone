CC = gcc
CFLAGS := -Wall -Wextra -Werror --pedantic --std=gnu11
LIBS := -lm

all: scone

scone: scone.c planets.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) $(LIBS) -o $@

planets.h: cfg.py kop.py rssk.cfg
	./kop.py > planets.h
