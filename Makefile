CC=gcc
CFLAGS=-O3 -Wall -Wextra -Werror -pedantic -Wconversion

binaries=euler runge_kutta

all: $(binaries)

euler: euler.c
	$(CC) $(CFLAGS) euler.c -o $@

runge_kutta: runge_kutta.c 
	$(CC) $(CFLAGS) runge_kutta.c -o $@

.PHONY: clean
clean:
	rm -f $(binaries) *.o *.dat
