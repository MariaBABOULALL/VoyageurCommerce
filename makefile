VC : Parseur.o
	gcc -o Parseur Parseur.o -lm
Parseur.o : Parseur.c
	gcc -c Parseur.c
