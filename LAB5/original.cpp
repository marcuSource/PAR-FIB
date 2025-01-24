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
#include <omp.h>

// huge pages
#include <sys/mman.h>
#include <sys/stat.h>

#include "fake.h"

// include files to visualize computed image
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

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
#define ROWS (5120)             /* size of display window in pixels */
#endif

#ifndef COLS
#define COLS (5120)             /* size of display window in pixels */
#endif

#ifndef TILE
#define TILE 32
#endif

#ifndef TILE_X
#define TILE_X 128
#endif

#define nTZOOMs 4
#define TZOOM (COLS/nTZOOMs)

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

void mandel_simple(int M[ROWS][COLS], double CminR, double CminI, double CmaxR, double CmaxI, double scale_real, double scale_imag, int maxiter)
{
      // Calcular
      for (int py = 0; py < ROWS; py++) {
	    for (int px = 0; px < COLS; px++) {
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
      }
}

void create_image(unsigned char *image, int matrix[ROWS][COLS], int height, int width, int maxiter)
{
      unsigned char filter[3][MAXITER];

      for (int i = 2; i < 3; i++) {
	    for (int j = 0; j < maxiter; j++)
		  filter[i][j] = rand();
	    for (int j = maxiter; j < MAXITER; j++)
		  filter[i][j] = 0;
      }
      for (int i = 0; i < 2; i++)
	    for (int j = 0; j < MAXITER; j++)
		  filter[i][j] = 0;

      for (int i = 0; i < 3; i++)
	    filter[i][maxiter] = 255;

      int ind = 0;
      int yy = height - 1;
      for (int y = 0; y < height; y++, yy--)
	    for (int x = 0; x < width; x++, ind += 3) {
		  ind = ((yy * width) + x) * 3;

		  image[ind] = filter[0][matrix[y][x]];
		  image[ind + 1] = filter[1][matrix[y][x]];
		  image[ind + 2] = filter[2][matrix[y][x]];
	    }
}

void create_grid(unsigned char *image, int height, int width)
{
      int ind;

      for (int i = 0; i < height; i++)
	    for (int j = TILE_X; j < width; j += TILE_X) {
		  ind = (i * width + j) * 3;
		  image[ind] = 255;
		  image[ind + 1] = 255;
		  image[ind + 2] = 255;
	    }

      for (int i = TILE; i < height; i += TILE)
	    for (int j = 0; j < width; j++) {
		  ind = (i * width + j) * 3;
		  image[ind] = 255;
		  image[ind + 1] = 255;
		  image[ind + 2] = 255;
	    }

      for (int i = 0; i < height; i++)
	    for (int j = TZOOM; j < width; j += TZOOM) {
		  ind = (i * width + j) * 3;
		  image[ind] = 255;
		  image[ind + 1] = 0;
		  image[ind + 2] = 0;
	    }

      for (int i = TZOOM; i < height; i += TZOOM)
	    for (int j = 0; j < width; j++) {
		  ind = (i * width + j) * 3;
		  image[ind] = 255;
		  image[ind + 1] = 0;
		  image[ind + 2] = 0;
	    }

}

#define HUGE_PAGE_SIZE (2 * 1024 * 1024)
#define ALIGN_TO_PAGE_SIZE(x) \
       (((x) + HUGE_PAGE_SIZE - 1) / HUGE_PAGE_SIZE * HUGE_PAGE_SIZE)

void *malloc_huge_pages(size_t size)
{
// Use 1 extra page to store allocation metadata
// (libhugetlbfs is more efficient in this regard)
      size_t real_size = ALIGN_TO_PAGE_SIZE(size + HUGE_PAGE_SIZE);
      char *ptr = (char *)mmap(NULL, real_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_HUGETLB, -1, 0);
      if (ptr == MAP_FAILED) {
	    // The mmap() call failed. Try to malloc instead
	    printf("malloc done\n");
	    ptr = (char *)malloc(real_size);
	    if (ptr == NULL)
		  return NULL;
	    real_size = 0;
      }
      // Save real_size since mmunmap() requires a size parameter
      *((size_t *)ptr) = real_size;
      // Skip the page with metadata
      return ptr + HUGE_PAGE_SIZE;
}

void free_huge_pages(void *ptr)
{
      if (ptr == NULL)
	    return;
      // Jump back to the page with metadata
      void *real_ptr = (char *)ptr - HUGE_PAGE_SIZE;

      // Read the original allocation size
      size_t real_size = *((size_t *)real_ptr);

      assert(real_size % HUGE_PAGE_SIZE == 0);
      if (real_size != 0)
            // The memory was allocated via mmap()
            // and must be deallocated via munmap()
	    munmap(real_ptr, real_size);
      else
	    // The memory was allocated via malloc()
	    // and must be deallocated via free()
	    free(real_ptr);
}

/*
*/

int main(int argc, char *argv[])
{
      int maxiter = 100;
      double real_min;
      double real_max;
      double imag_min;
      double imag_max;
      int width = COLS;         // dimensions of display window
      int height = ROWS;
      double size = N, x0 = 0, y0 = 0;
      int *Hmatrix;
      unsigned char *image;

      TRACE_START_MARK;
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
		  fprintf(stderr, "       -o to write computed image (mandel_image.jpg) and histogram (mandel_histogram.out) to disk (default no file generated)\n");
		  fprintf(stderr, "       -h to produce histogram of values in computed image (default no histogream)\n");
		  fprintf(stderr, "       -d to display computed image (default no display)\n");
		  fprintf(stderr, "       -i to specify maximum number of iterations at each point (default 100)\n");
		  fprintf(stderr, "       -c to specify the center x0+iy0 of the square to compute (default origin)\n");
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

      if (output2file) {
	    image = (unsigned char *)malloc(width * height * 3);
	    if (image == NULL) {
		  printf("Malloc (image) ERROR\n");
		  exit(0);
	    }
      }

      Hmatrix = (int *)malloc_huge_pages(height * width * sizeof(int));
      if (Hmatrix == NULL) {
	    printf("Malloc (Hmatrix) ERROR\n");
	    exit(0);
      }

      if (output2histogram)
	    histogram = (int *)calloc(maxiter, sizeof(int));

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

      mandel_simple((int (*)[COLS])Hmatrix, real_min, imag_min, real_max, imag_max, scale_real, scale_imag, maxiter);

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
	    create_image(image, (int (*)[COLS])Hmatrix, height, width, maxiter);
	    //create_grid (image, height, width);
	    stbi_write_png("mandel_image.jpg", width, height, 3, image, 0);
	    if (output2histogram)
		  if (fwrite(histogram, sizeof(int), maxiter, fp) != maxiter)
			fprintf(stderr, "Error when writing histogram to file\n");
      }

      free_huge_pages(Hmatrix);
      Hmatrix = NULL;
      if (output2file) {
	    free(image);
	    image = NULL;
      }

      if (output2display) {
	    /* Wait for user response, then exit program */
	    if (setup_return == EXIT_SUCCESS) {
		  interact(display, &win, width, height, real_min, real_max, imag_min, imag_max);
	    }
      }
      TRACE_END_MARK;
      return EXIT_SUCCESS;
}
