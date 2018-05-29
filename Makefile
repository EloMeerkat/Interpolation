CC=gcc
CFLAGS=-I. -lm
DEPS = methodes.h fonctions.h
OBJ = main.o methodes.o fonctions.o

%.o:%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)
