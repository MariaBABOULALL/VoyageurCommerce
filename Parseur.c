/*
 * Parseur.c
 *
 *  Created on: 14 mars 2014
 *      Author: ymoreno
 */

#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>

#define MAXMOT 256
#define MAXS 500

/**
 * Structure pour representer un cycle
 */
typedef struct t_cycle
{
  int taille;   //la taille du cycle en construction
  double poids; //le co�t du cycle
  int c[MAXS];  //liste des "taille" sommets
} t_cycle;


/**
 * Charge le CSV des coordonn�es des villes.
 *
 * @param [in] f le fichier
 * @param [out] nb_villes le nombre de villes de l'instance
 * @param [out] dist le tableau des nb_villes*nb_villes distances
 * @param [out] absc le tableau des abscisses des villes
 * @param [out] ord le tableau des ordonn�es des villes
 */
void lire_donnees(const char *f, unsigned int *nb_villes, double ***dist,  double **absc, double **ord)
{
  //double *absc; ///tableau des ordonn�es
  //double *ord;  /// tableau des abscisses
  char ligne[MAXMOT];
  FILE * fin = fopen(f,"r");

  if(fin != NULL)
    {
      //On recupere le nombre de villes
      fgets(ligne, MAXMOT, fin);
      *nb_villes = atoi(ligne);
      (*dist) = (double**)malloc(*nb_villes * sizeof(double*));
      (*absc) = (double*)malloc(*nb_villes * sizeof(double));
      (*ord) = (double*)malloc(*nb_villes * sizeof(double));
      int i = 0;
      while (fgets(ligne, MAXMOT, fin) != NULL)
	{
	  char *p = strchr(ligne, ';');
	  ligne[strlen(ligne) - strlen(p)]='\0';
	  p = &p[1];
	  (*absc)[i] = atof(ligne);
	  (*ord)[i] = atof(p);
	  i = i + 1;
	}
    }
  else
    {
      printf("Erreur de lecture du fichier.\n");
      exit(2);
    }
  fclose(fin);
  int i,j;

  //Calclul des distances
  for(i = 0; i < *nb_villes; i++)
    {
      (*dist)[i] = (double*) malloc(*nb_villes * sizeof(double));
      for(j = 0; j < *nb_villes; j++)
	{
	  (*dist)[i][j] = sqrt( ((*absc)[i] - (*absc)[j])* ((*absc)[i] - (*absc)[j]) + ((*ord)[i] - (*ord)[j]) * ((*ord)[i] - (*ord)[j]) );
	}
    }
}


/**
 * Supprime la structure des distances
 *
 * @param [in] nb_villes le nombre de villes.
 * @param [in,out] distances le tableau � supprimer.
 * @param [in,out] abscisses un autre tableau � supprimer.
 * @param [in,out] ordonnees encore un autre tableau � supprimer.
 */
void supprimer_distances_et_coordonnees(const int nb_villes, double **distances, double *abscisses, double *ordonnees)
{
  int i;
  for(i = 0; i < nb_villes; i++)
    {
      free(distances[i]);
    }
 free(distances);
 free(abscisses);
 free(ordonnees);
}


/**
 * Export le cycle dans un fichier HTML pour pouvoir �tre visualis�
 * dans l'applet.
 *
 * @param [in] cycle le cycle � afficher
 */
void afficher_cycle_html(const t_cycle cycle, double *posX, double *posY)
{
  FILE * fout = fopen("DisplayTsp2.html","w");
  if(fout != NULL)
    {
      int i;
      fprintf(fout, "<html>\n <applet codebase=\".\" code=\"DisplayTsp.class\" width=200 height=200>\n");
      fprintf(fout, "<param name = Problem value = \"custom\">\n");
      fprintf(fout, "<param name = Problem CitiesPosX value = \"");
      for(i = 0; i < cycle.taille; i++)
	fprintf(fout,"%f;",posX[i]);
      fprintf(fout, "\">\n");
      fprintf(fout, "<param name = Problem CitiesPosY value = \"");
      for(i = 0; i < cycle.taille; i++)
	fprintf(fout,"%f;",posY[i]);
      fprintf(fout, "\">\n");
      fprintf(fout, "<param name = Parcours value = \"");
      fprintf(fout,"%d",cycle.c[0]);
      for(i = 1; i < cycle.taille; i++)
	fprintf(fout,"-%d",cycle.c[i]);
      fprintf(fout,"\">\n</applet>\n </html>\n");
    }
  fclose(fout);
}

/**
 * Affiche le tableau des distances.
 *
 * @param [in] nb le nombre de villes
 * @param [in] distances le tableau
 */
void afficher_distances(const int nb, double **distances)
{
  unsigned int i ;
  unsigned int j ;
  for(i = 0  ; i < nb; i++)
    {
      for(j = 0 ; j < nb ; j++)
	printf("%f ", distances[i][j]);
      printf("\n");
    }
  printf("\n");
}


/**
 * Fonction de comparaison pour le trie des ar�tes par leur poids.
 *
 * @param [in] v1 pointeur vers un triplet (i,j,poids)
 * @param [in] v2 pointeur verts un triplet (i,j,poids)
 * @return vrai si poid v1 < poids v2
 */
int comparer(const void *v1, const void *v2)
{
  double **px1 = (double **) v1;
  double **px2 = (double **) v2;

  double *x1 = *px1;
  double *x2 = *px2;
  if(x1[2] - x2[2] < 0)
    return -1;
  else
    {
      if(x1[2] - x2[2] == 0)
	return 0;
      else
	return 1;
    }
}

/**
 * Construit un tableau de n*(n-1)/2 ar�tes tri� selon leur poids.
 *
 * @note utile pour le Kruskal
 *
 * @param [in] n le nombre de villes
 * @param [in] d tableau des n x n distances.
 * @return tableau d'ar�tes tri�es T[i][j] = poids_ij
 */
double **trier_aretes(const int n, double **d)
{
  assert(d);

  int nb_aretes = n * (n - 1) / 2;
  double **T = (double **)malloc(nb_aretes * sizeof(double *));
  int i, j;
  int a = 0;

  //On initialise la structure d'ar�tes
  for(i = 0; i < n-1; i++)
    {
      for(j = i+1; j < n; j++)
	{
	  T[a] = (double *)malloc(3 * sizeof(double));
	  T[a][0] = i;
	  T[a][1] = j;
	  T[a][2] = d[i][j];
	  a++;
	}
    }

  //Appel au quicksort avec la bonne fonction de comparaison
  qsort(T, a, sizeof(T[0]), comparer);


  //Decommenter pour v�rifier le tri
  /*
  for(i = 0; i < a; i++)
    printf("%f ", T[i][2]);
  printf("\n");
  */
  return T;
}


/**
 * Supprime le tableau des ar�tes.
 *
 * @param [in] nb_villes le nombre de villes
 * @param [in,out] T le tableau � supprimer
 */
void supprimer_aretes(const int nb_villes, double **T)
{
  assert(T);

  int nb_aretes = nb_villes*(nb_villes - 1 ) / 2;
  unsigned int i;
  for( i = 0; i < nb_aretes ; ++i)
    free( T[i] );

  free(T);
}

/**
 * Fonction de comparaison de villes et chemins
 *
 * @param [ville] ville à vérifier
 * @param [chemin_courant] pointeur vers le chemin à verifier
 * @return faux(1) si ville est présente dans chemin
 */
int est_dans_chemin(int ville,t_cycle *chemin_courant)
{
	int resultat = 0;
  for (int k=0;k<chemin_courant->taille;k++)
  {
	  if (ville == chemin_courant->c[k])
	  {
		  resultat = 1;
		  return resultat;
	  }
  }
  return resultat;
}

/**
 * Fonction de remplacement de deux chemins (cycles)
 *
 * @param [chemin_a_remplacer]
 * @param [chemin_a_copier]
 */
void remplacement_tcycles(t_cycle* chemin_a_remplacer, t_cycle* chemin_a_copier)
{
	chemin_a_remplacer->poids=chemin_a_copier->poids;
	chemin_a_remplacer->taille=chemin_a_copier->taille;
	for (int k=0;k<chemin_a_copier->taille;k++)
	{
		chemin_a_remplacer->c[k]=chemin_a_copier->c[k];
	}
}

/**
 * Fonction d'ajout de ville
 *
 * @param [ville] ville a ajouter
 * @param [chemin] chemin de destination
 * @param [dist] tableau de distances
 */
void ajouter_ville(int ville, t_cycle *chemin, double *** dist)
{
	 // distance entre la dernière ville du chemin et la prochaine ville a ajouter
	chemin->poids+=dist[chemin->c[(chemin->taille)-1]][ville];
	 // augmenter la taille du tableau puis rajouter la ville
	chemin->c[(chemin->taille)]=ville;
	chemin->taille++;
}

/**
 * Fonction de retrait
 *
 * @param [ville] ville a retirer
 * @param [chemin] chemin de destination
 * @param [dist] tableau de distances
 */
void retirer_ville(int ville, t_cycle *chemin, double *** dist)
{
	 // distance entre la dernière ville du chemin et la prochaine ville a retirer
	chemin->poids-=dist[chemin->c[(chemin->taille)-1]][ville];
	 // reduire la taille du tableau puis retirer la ville
	chemin->c[(chemin->taille)]=0;
	chemin->taille--;
}

void PVC_EXACT_NAIF(int nb_villes, double *** dist, t_cycle* chemin, t_cycle* meilleur)
{
	if (chemin->taille == nb_villes)
	{
		if (chemin->poids < meilleur->poids)
		{
			remplacement_tcycles(meilleur,chemin);
		}
	}
	else
	{
		for (int ville = 0; ville<nb_villes; ville ++)
		{
			ajouter_ville(ville,chemin,dist);
			PVC_EXACT_NAIF(nb_villes,dist,chemin,meilleur);
			retirer_ville(ville,chemin,dist);
		}
	}
}

/**
 * Fonction main.
 */
int main (int argc, char *argv[])
{
  double **distances;
  double *abscisses;
  double *ordonnees;
  unsigned int nb_villes;

  //Initialisation du timer pour mesurer des temps (compiler avec -lrt)
  struct timespec myTimerStart;
  clock_gettime(CLOCK_REALTIME, &myTimerStart);

  //Exemple de mesure du temps
  lire_donnees("defi250.csv", &nb_villes, &distances, &abscisses, &ordonnees);

  //R�cup�ration du timer et affichage
  struct timespec current;
  clock_gettime(CLOCK_REALTIME, &current); //Linux gettime
  double elapsed_in_ms =    (( current.tv_sec - myTimerStart.tv_sec) *1000 +
          ( current.tv_nsec - myTimerStart.tv_nsec)/1000000.0);
  printf("Temps pass� (ms) : %lf\n", elapsed_in_ms);


  //Affichage des distances
  //afficher_distances(nb_villes,distances);

  //naif
  t_cycle cycle;
  cycle.taille=3;
  cycle.c[0]=0;
  cycle.c[1]=1;
  cycle.c[2]=2;
  afficher_cycle_html(cycle, abscisses, ordonnees);

  double ** Aretes =  trier_aretes(nb_villes, distances);
  /// <-- Kruskal Here
  supprimer_aretes(nb_villes, Aretes);

  supprimer_distances_et_coordonnees(nb_villes, distances, abscisses, ordonnees);
  return 0;
}



