#CFLAGS regroupe les options de compilation (Généralement cette variable est nommées CFLAGS pour une compilation en C, CXXFLAGS pour le C++).
#LDFLAGS regroupe les options de l'édition de liens.
#EXEC contient le nom des exécutables à générer.


CC=gcc
CFLAGS=-Wall
LDFLAGS=-lm -lgsl -lgslcblas
EXEC=a.out
SRC= $(wildcard *.c)
OBJ= $(SRC:.c=.o)


$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

$(OBJ): $(wildcard *.h)
#on fait ca car le .h a des flag et donc si on modifie le .h il faut tt recompiler.

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
