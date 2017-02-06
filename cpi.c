/* ------------------------------
   $Id: mandel-seq.c,v 1.2 2008/03/04 09:52:55 marquet Exp $
   ------------------------------------------------------------

   Affichage de l'ensemble de Mandelbrot.
   Version sequentielle.

*/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "string.h"

/* Valeur par defaut des parametres */
#define N_ITER  255		/* nombre d'iterations */

#define X_MIN   -1.78		/* ensemble de Mandelbrot */
#define X_MAX   0.78
#define Y_MIN   -0.961
#define Y_MAX   0.961

#define X_SIZE  4096		/* dimension image */
#define Y_SIZE  3072
#define FILENAME "mandel.ppm"	/* image resultat */

typedef struct {
    int x_size, y_size;		/* dimensions */
    char *pixels;		/* matrice linearisee de pixels */
} picture_t;

static void
usage()
{
    fprintf(stderr, "usage : ./mandel [options]\n\n");
    fprintf(stderr, "Options \t Signification \t\t Val. defaut\n\n");
    fprintf(stderr, "-n \t\t Nbre iter. \t\t %d\n", N_ITER);
    fprintf(stderr, "-b \t\t Bornes \t\t %f %f %f %f\n",
	    X_MIN, X_MAX, Y_MIN, Y_MAX);
    fprintf(stderr, "-d \t\t Dimensions \t\t %d %d\n", X_SIZE, Y_SIZE);
    fprintf(stderr, "-f \t\t Fichier \t\t %s\n", FILENAME);

    exit(EXIT_FAILURE);
}

static void
parse_argv (int argc, char *argv[],
	    int *n_iter,
	    double *x_min, double *x_max, double *y_min, double *y_max,
	    int *x_size, int *y_size,
	    char **path)
{
    const char *opt = "b:d:n:f:";
    int c;

    /* Valeurs par defaut */
    *n_iter = N_ITER;
    *x_min  = X_MIN;
    *x_max  = X_MAX;
    *y_min  = Y_MIN;
    *y_max  = Y_MAX;
    *x_size = X_SIZE;
    *y_size = Y_SIZE;
    *path   = FILENAME;

    /* Analyse arguments */
    while ((c = getopt(argc, argv, opt)) != EOF) {
	switch (c) {
	    case 'b': 		/* domaine */
		sscanf(optarg, "%lf", x_min);
		sscanf(argv[optind++], "%lf", x_max);
		sscanf(argv[optind++], "%lf", y_min);
		sscanf(argv[optind++], "%lf", y_max);
		break;
	    case 'd':		/* largeur hauteur */
		sscanf(optarg, "%d", x_size);
		sscanf(argv[optind++], "%d", y_size);
		break;
	    case 'n':		/* nombre d'iterations */
		*n_iter = atoi(optarg);
		break;
	    case 'f':		/* fichier de sortie */
		*path = optarg;
		break;
	    default:
		usage();
	}
    }
}

static void
init_picture (picture_t *pict, int x_size, int y_size)
{
    pict->y_size = y_size;
    pict->x_size = x_size;
    pict->pixels = malloc(y_size * x_size); /* allocation espace memoire */
}

/* Enregistrement de l'image au format ASCII .ppm */
static void
save_picture (const picture_t *pict, const char *pathname)
{
    unsigned i;
    FILE *f = fopen(pathname, "w");

    fprintf(f, "P6\n%d %d\n255\n", pict->x_size, pict->y_size);
    for (i = 0 ; i < pict->x_size * pict->y_size; i++) {
	char c = pict->pixels[i];
	fprintf(f, "%c%c%c", c, c, c); /* monochrome blanc */
    }

    fclose (f);
}

static void
compute (picture_t *pict,
	 int nb_iter,
	 double x_min, double x_max, double y_min, double y_max, int myid, int numprocs)
{
    int pos = 0;
    int iy, ix, i;
    double pasx = (x_max - x_min) / pict->x_size, /* discretisation */
	   pasy = (y_max - y_min) / pict->y_size;
     printf("max = %f, min =%f\n", y_max, y_min);
     printf("pas = %f, max-min = %f\n", pasy, (y_max - y_min));
     printf(">myid = %d, pic_dimensions (%d x %d)\n", myid, pict->x_size, pict->y_size);

    /* Calcul en chaque point de l'image */
    for (iy = 0 ; iy < pict->y_size ; iy++) {
	     for (ix = 0 ; ix < pict->x_size; ix++) {
	        double a = x_min + ix * pasx, b = y_max - iy * pasy, x = 0, y = 0;
	        for (i = 0 ; i < nb_iter ; i++) {
		          double tmp = x;
		          x = x * x - y * y + a;
		          y = 2 * tmp * y + b;
		          if (x * x + y * y > 4) /* divergence ! */
		            break;
	        }
       pict->pixels[pos++] = (double) i / nb_iter * 255;
	    }
    }
    printf(">myid = %d and just finish compute\n", myid);
}

int
main (int argc, char *argv[])
{

    int n_iter,			/* degre de nettete  */
	x_size, y_size;		/* & dimensions de l'image */
    double x_min, x_max, y_min, y_max; /* bornes de la representation */
    char *pathname;		/* fichier destination */


	int numprocs, myid, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];


    parse_argv(argc, argv,
	       &n_iter,
	       &x_min, &x_max, &y_min, &y_max,
	       &x_size, &y_size, &pathname);


         MPI_Init(&argc,&argv);
           MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
           MPI_Comm_rank(MPI_COMM_WORLD,&myid);
           MPI_Get_processor_name(processor_name,&namelen);
       picture_t pict;

         double part_proc = (y_max - y_min) / numprocs;

         printf("I am task %d out of %d\n" ,myid, numprocs);
         printf("> myid = %d part_proc = %f\n", myid, part_proc);

    // init_picture (& pict, x_size, y_size);
    // compute (& pict, n_iter, x_min, x_max, y_min, y_max);
    init_picture (& pict, x_size, y_size/numprocs);

    printf("> myid = %d picture initialized\n", myid);
    printf("> pic addr = %p\n", &pict);

    compute (& pict, n_iter, x_min, x_max, y_min + ((numprocs - myid - 1) * part_proc), (y_max - (myid * part_proc)), myid, numprocs);

    printf("> I am process %d and I finished my compute, I start saving my pic\n", myid);


    char str[15];
    sprintf(str, "%d", myid);
    strcat(&pathname, &str);

    save_picture (&pict, &pathname);

    MPI_Finalize();

    exit(EXIT_SUCCESS);
}
