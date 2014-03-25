VC : Parseur.o
	gcc -o Parseur Parseur.o -lm -lrt
Parseur.o : Parseur.c
	gcc -c Parseur.c
