/*
 * Sequential Mandelbrot program
 *
 * This program computes and displays all or part of the Mandelbrot
 * set. By default, it examines all points in the complex plane
 * that have both real and imaginary parts between -2 and 2.
 * Command-line parameters allow zooming in on a specific part of
 * this range.
 *
 * Basic usage:
 *   mandel [-d -i maxiter -c x0 y0 -s size -w windowsize]
 * where
 *   maxiter denotes the maximum number of iterations at each point -- by default 100
 *   x0, y0, and size specify the range to examine (a square
 *     centered at (x0 + iy0) of size 2*size by 2*size -- by default,
 *     a square of size 4 by 4 centered at the origin)
 *   windowsize denotes the size of the image (diplay window) to compute
 *
 * Output: if -d option is used, a graphical display as described in
 *         Wilkinson & Allen, is displayed using the X Window system,
 *         plus text output to standard output showing the above parameters,
 *         plus execution time in seconds.
 *
 * Additional parameters:
 *   If -h option is used, the program computes the histogram of values in the image
 *   If -o option is used, the program saves output and histogram (if computed) to file
 *
 * Code based on the original code from Web site for Wilkinson and Allen's
 * text on parallel programming:
 * http://www.cs.uncc.edu/~abw/parallel/par_prog/
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <malloc.h>

// include file to compile with Tareador
#include <tareador.h>

// include files to visualize computed image
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

// include file and stuff to measure execution times
#include <sys/time.h>
double getusec_()
{
      struct timeval time;
      gettimeofday(&time, NULL);
      return ((double)time.tv_sec * (double)1e6 + (double)time.tv_usec);
}

double stamp;
#define START_COUNT_TIME stamp = getusec_();
#define STOP_COUNT_TIME(_m) stamp = getusec_() - stamp;\
                        stamp = stamp/1e6;\
                        printf ("%s %0.6f\n",(_m), stamp);

// Default values for things
#define N           2           /* size of problem space (x, y from -N to N) */
#ifndef ROWS
#define ROWS 32                 /* size of display window in pixels */
#endif

#ifndef COLS
#define COLS 32                 /* size of display window in pixels */
#endif

#ifndef TILE
#define TILE 8
#endif

// Structure definition for complex numbers
typedef struct {
      double real, imag;
} complex;

// Functions for GUI
#include "mandelbrot-gui.h"     /* has setup(), interact() */

// Global variables to output results
// output to file
int output2file = 0;
FILE *fp = NULL;

// output to display
int output2display = 0;
Display *display;
Window win;
GC gc;
int setup_return;
long min_color = 0, max_color = 0;
double scale_color;
double scale_real, scale_imag;

// output as histogram
int output2histogram = 0;
int *histogram;

// dummy parameter
int user_param = 1;

#define MAXITER 10000

int pixel_dwell(int width, int height, double CminR, double CminI, double CmaxR, double CmaxI, int x, int y, double scale_real, double scale_imag, int maxiter)
{
      int dwell = 0;

      double cR, cI, zR, zI;

      cR = CminR + ((double)x * scale_real);
      cI = CminI + (double)(height - 1 - y) * scale_imag;

      // z = c;
      zR = cR;
      zI = cI;

      double tmpR, tmpI;
      double lengthsq;
      do {
	    //z = z * z + c;
	    tmpR = zR * zR - zI * zI + cR;
	    tmpI = 2 * zR * zI + cI;
	    zR = tmpR;
	    zI = tmpI;
	    lengthsq = (zR * zR + zI * zI);
	    dwell++;
      } while ((dwell < maxiter) && (lengthsq < (N * N)));

      return dwell;
}

void mandel_tiled_rec(int M[ROWS][COLS], int NRows, int NCols, int start_fil, int start_col, double CminR, double CminI, double CmaxR, double CmaxI, double scale_real, double scale_imag, int maxiter)
{
      int non_first_call= (NRows!=ROWS) | (NCols!=COLS);
      int equal;
      int x, y;

      y = start_fil;
      x = start_col;
      equal = non_first_call;
      //check horizontal borders
      for (int px = x; px < x + NCols; px++) {
	    M[y][px] = pixel_dwell(COLS, ROWS, CminR, CminI, CmaxR, CmaxI, px, y, scale_real, scale_imag, maxiter);
	    M[y + NRows - 1][px] = pixel_dwell(COLS, ROWS, CminR, CminI, CmaxR, CmaxI, px, y + NRows - 1, scale_real, scale_imag, maxiter);
	    equal = equal & (M[y][x] == M[y][px]);
	    equal = equal & (M[y][x] == M[y + NRows - 1][px]);
      }
      //check vertical borders
      for (int py = y; py < y + NRows; py++) {
	    M[py][x] = pixel_dwell(COLS, ROWS, CminR, CminI, CmaxR, CmaxI, x, py, scale_real, scale_imag, maxiter);
	    M[py][x + NCols - 1] = pixel_dwell(COLS, ROWS, CminR, CminI, CmaxR, CmaxI, x + NCols - 1, py, scale_real, scale_imag, maxiter);
	    equal = equal & (M[y][x] == M[py][x]);
	    equal = equal & (M[y][x] == M[py][x + NCols - 1]);
      }
      //check if we can avoid computation of tile
      if (equal) {
	    // fill all with the same value: base case
		long color = (long)((M[y][x] - 1) * scale_color) + min_color;
		if (output2histogram)
		   histogram[M[y][x] - 1]+= (NRows*NCols);
	    for (int py = y; py < y + NRows; py++)
		  for (int px = x; px < x + NCols; px++) {
			M[py][px] = M[y][x];
			if (output2display) {
			      /* Scale color and display point  */
			      if (setup_return == EXIT_SUCCESS) {
				    XSetForeground(display, gc, color);
				    XDrawPoint(display, win, gc, px, py);
			      }
			}
		  }
      } else {
	    if (NCols <= TILE) {

		tareador_start_task("leaf");
		  // computation of tile: base case
		  for (int py = y; py < y + NRows; py++)
			for (int px = x; px < x + NCols; px++) {
			      M[py][px] = pixel_dwell(COLS, ROWS, CminR, CminI, CmaxR, CmaxI, px, py, scale_real, scale_imag, maxiter);
			      if (output2histogram)
				    histogram[M[py][px] - 1]++;
			      if (output2display) {
				    /* Scale color and display point  */
				    long color = (long)((M[py][px] - 1) * scale_color) + min_color;
				    if (setup_return == EXIT_SUCCESS) {
					  XSetForeground(display, gc, color);
					  XDrawPoint(display, win, gc, px, py);
				    }
			      }
			}
		  tareador_end_task("leaf");
	    } else {
		  // computation of tile: recursive cases
		  if (NRows > TILE) {
			mandel_tiled_rec(M, NRows / 2, NCols / 2, start_fil, start_col, CminR, CminI, CmaxR, CmaxI, scale_real, scale_imag, maxiter);
			mandel_tiled_rec(M, NRows / 2, NCols / 2, start_fil, start_col + NCols / 2, CminR, CminI, CmaxR, CmaxI, scale_real, scale_imag, maxiter);
			mandel_tiled_rec(M, NRows / 2, NCols / 2, start_fil + NRows / 2, start_col, CminR, CminI, CmaxR, CmaxI, scale_real, scale_imag, maxiter);
			mandel_tiled_rec(M, NRows / 2, NCols / 2, start_fil + NRows / 2, start_col + NCols / 2, CminR, CminI, CmaxR, CmaxI, scale_real, scale_imag, maxiter);
		  } else {
			mandel_tiled_rec(M, NRows, NCols / 2, start_fil, start_col, CminR, CminI, CmaxR, CmaxI, scale_real, scale_imag, maxiter);
			mandel_tiled_rec(M, NRows, NCols / 2, start_fil, start_col + NCols / 2, CminR, CminI, CmaxR, CmaxI, scale_real, scale_imag, maxiter);
		  }

	    }
      }
}

int main(int argc, char *argv[])
{
      int maxiter = 100;
      double real_min;
      double real_max;
      double imag_min;
      double imag_max;
      int width = COLS;         // dimensions of display window
      int height = ROWS;
      double size = N, x0 = -1.81641, y0= -0.203125;
      int *Hmatrix;

      // Process command-line arguments
      for (int i = 1; i < argc; i++) {
	    if (strcmp(argv[i], "-d") == 0) {
		  output2display = 1;
	    } else if (strcmp(argv[i], "-h") == 0) {
		  output2histogram = 1;
	    } else if (strcmp(argv[i], "-i") == 0) {
		  maxiter = atoi(argv[++i]);
	    } else if (strcmp(argv[i], "-c") == 0) {
		  x0 = atof(argv[++i]);
		  y0 = atof(argv[++i]);
	    } else if (strcmp(argv[i], "-s") == 0) {
		  size = atof(argv[++i]);
	    } else if (strcmp(argv[i], "-u") == 0) {
		  user_param = atof(argv[++i]);
	    } else if (strcmp(argv[i], "-o") == 0) {
		  output2file = 1;
		  if ((fp = fopen("mandel_histogram.out", "wb")) == NULL) {
			fprintf(stderr, "Unable to open file\n");
			return EXIT_FAILURE;
		  }
	    } else {
		  fprintf(stderr, "Usage: %s [-o -h -d -i maxiter -c x0 y0 -s size]\n", argv[0]);
		  fprintf(stderr, "       -o to write the computed histogram (mandel_histogram.out) to disk (default no file generated)\n");
		  fprintf(stderr, "       -h to produce histogram of values in computed image (default no histogream)\n");
		  fprintf(stderr, "       -d to display computed image (default no display)\n");
		  fprintf(stderr, "       -i to specify maximum number of iterations at each point (default 100)\n");
		  fprintf(stderr, "       -c to specify the center x0+iy0 of the square to compute\n");
		  fprintf(stderr, "       -s to specify the size of the square to compute (default 2, i.e. size 4 by 4)\n");
		  return EXIT_FAILURE;
	    }
      }
      real_min = x0 - size;
      real_max = x0 + size;
      imag_min = y0 - size;
      imag_max = y0 + size;
      width = COLS;
      height = ROWS;

      // Produce text output
      fprintf(stdout, "\n");
      fprintf(stdout, "Computation of the Mandelbrot set with:\n");
      fprintf(stdout, "    center = (%g, %g) \n    size = %g\n    maximum iterations = %d\n", (real_max + real_min) / 2, (imag_max + imag_min) / 2, (real_max - real_min) / 2, maxiter);
      fprintf(stdout, "\n");

      Hmatrix = (int *)malloc(height * width * sizeof(int));
      if (Hmatrix == NULL) {
	    printf("Malloc (Hmatrix) ERROR\n");
	    exit(0);
      }

      if (output2histogram)
	    histogram = calloc(maxiter, sizeof(int));

      if (output2display) {
	    /* Initialize for graphical display */
	    setup_return = setup(width, height, &display, &win, &gc, &min_color, &max_color);
	    if (setup_return != EXIT_SUCCESS) {
		  fprintf(stderr, "Unable to initialize display, continuing\n");
		  return EXIT_FAILURE;
	    }
      }

      // Compute factors to scale computational region to window
      scale_real = (double)(real_max - real_min) / (double)width;
      scale_imag = (double)(imag_max - imag_min) / (double)height;

      if (output2display) {
	    // Compute factor for color scaling
	    scale_color = (double)(max_color - min_color) / (double)(maxiter - 1);
      }

      // Start timing
      if (!output2display) {
	    START_COUNT_TIME;
      }

      tareador_ON();
      mandel_tiled_rec((int (*)[width])Hmatrix, height, width, 0, 0, real_min, imag_min, real_max, imag_max, scale_real, scale_imag, maxiter);
      tareador_OFF();

      // End timing
      if (!output2display) {
	    STOP_COUNT_TIME("Total execution time (in seconds):");
	    fprintf(stdout, "\n");
      }

      fprintf(stdout, "Mandelbrot set: Computed\n");
      if (output2histogram)
	    fprintf(stdout, "Histogram for Mandelbrot set: Computed\n");
      else
	    fprintf(stdout, "Histogram for Mandelbrot set: Not computed\n");

      // Make sure all output is written
      if ((output2display) && (setup_return == EXIT_SUCCESS))
	    XFlush(display);

      if ((output2file) && (fp != NULL)) {
	    fprintf(stdout, "Writing output file to disk: mandel_image.jpg\n");
	    if (output2histogram)
		  if (fwrite(histogram, sizeof(int), maxiter, fp) != maxiter)
			fprintf(stderr, "Error when writing histogram to file\n");
      }

      free(Hmatrix);
      Hmatrix = NULL;

      if (output2display) {
	    /* Wait for user response, then exit program */
	    if (setup_return == EXIT_SUCCESS) {
		  interact(display, &win, width, height, real_min, real_max, imag_min, imag_max);
	    }
	    return EXIT_SUCCESS;
      }
}
