#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char * argv[])
{
  FILE * fp = fopen ("conf.gro", "w");
  double d = 2;
  int nx, ny, nz;
  nx = ny = nz = 6;

  fprintf (fp, "\n%d\n", nx*ny*nz);

  int count = 1;
  for (unsigned i = 0; i < nx; ++i){
    for (unsigned j = 0; j < ny; ++j){
      for (unsigned k = 0; k < nz; ++k){
	fprintf(fp, "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
		count, "ljr", "lja", count,
		i * d + 0.5 * d ,
		j * d + 0.5 * d + i * 0.16,
		k * d + 0.5 * d + j * 0.16 * i * 0.16,
		0., 0., 0.);
      }
    }
  }

  fprintf (fp, "%f %f %f\n", nx*d, ny*d, nz*d);
  
  fclose(fp);
}
