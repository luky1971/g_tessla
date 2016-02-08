CC = cc

# add -DINEXACT=volatile to CFLAGS if your processor uses internal extended precision floating-point registers

CFLAGS = -O3

predicates.o: predicates.c
	$(CC) $(CFLAGS) -o predicates.o -c predicates.c

clean:
	rm predicates.o
