#define CPP_FILE
#include "GromacsFileManager.h"
#include <stdio.h>
#include <string.h>

void
GromacsFileManager::readGroFile (const char * filename,
				 IndexType * resdindex,
				 char * resdname,
				 char * atomname,
				 IndexType * atomindex,
				 ScalorType * posix,
				 ScalorType * posiy,
				 ScalorType * posiz,
				 ScalorType * velox,
				 ScalorType * veloy,
				 ScalorType * veloz,
				 ScalorType * boxx,  
				 ScalorType * boxy,  
				 ScalorType * boxz)
{
  FILE * fp = fopen (filename, "r");
  if (fp == NULL){
    throw MDExcptCannotOpenFile (filename);
  }
  char line[1024];
  fgets (line, 1024, fp);
  IndexType tmpn;
  fgets (line, 1024, fp);
  if (sscanf (line, "%d", &tmpn) != 1){
    throw MDExcptWrongFileFormat ("GromacsFileManager::readGroFile", filename);
  }

  int n1, n2;
  char s1[10], s2[10];
  ScalorType a, b, c;
  ScalorType d, e, f;
  char tmp[8];
  
  for (IndexType i = 0; i < tmpn; ++i){
    fgets (line, 1024, fp);
    for (IndexType j = 0; j < 5; ++j){
      tmp[j] = line[j+5];
    }
    tmp[5] = '\0';
    if (sscanf (tmp, "%s", s1) != 1){
      throw MDExcptWrongFileFormat ("GromacsFileManager::readGroFile", filename);
    }
    for (IndexType j = 0; j < 5; ++j){
      tmp[j] = line[j+10];
    }
    tmp[5] = '\0';
    if (sscanf (tmp, "%s", s2) != 1){
      throw MDExcptWrongFileFormat ("GromacsFileManager::readGroFile", filename);
    }
    int tag = sscanf (line, "%5d", &n1);
    if (tag != 1){
      throw MDExcptWrongFileFormat ("GromacsFileManager::readGroFile", filename);
    }
    tag = sscanf (&line[15], "%5d%f%f%f%f%f%f", &n2, &a, &b, &c, &d, &e, &f);
    switch (tag){
    case 7:
	velox[i] = d;
	veloy[i] = e;
	veloz[i] = f;
	break;
    case 4:
	velox[i] = 0.f;
	veloy[i] = 0.f;
	veloz[i] = 0.f;
	break;
    default:
	throw MDExcptWrongFileFormat ("GromacsFileManager::readGroFile", filename);
    }
    strcpy (&resdname[i*StringSize], s1);
    strcpy (&atomname[i*StringSize], s2);
    resdindex[i] = n1;
    atomindex[i] = n2;
    posix[i] = a;
    posiy[i] = b;
    posiz[i] = c;
  }
  fgets (line, 1024, fp);
  int tag = sscanf (line, "%f%f%f", &a, &b, &c);
  if (tag != 3){
    throw MDExcptWrongFileFormat ("GromacsFileManager::readGroFile", filename);
  }
  *boxx = a;
  *boxy = b;
  *boxz = c;
  
  fclose (fp);
}


void
GromacsFileManager::writeGroFile (FILE * fp,
				  const IndexType num,
				  const IndexType * resdindex,
				  const char * resdname,
				  const char * atomname,
				  const IndexType * atomindex,
				  const ScalorType * posix,
				  const ScalorType * posiy,
				  const ScalorType * posiz,
				  const ScalorType * velox,
				  const ScalorType * veloy,
				  const ScalorType * veloz,
				  const ScalorType boxx,  
				  const ScalorType boxy,  
				  const ScalorType boxz)
{
  fprintf (fp, "\n%d\n", num);
  for (IndexType i = 0; i < num; ++i){
    fprintf(fp, "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
	    resdindex[i], &resdname[i*StringSize],
	    &atomname[i*StringSize], atomindex[i],
	    posix[i], posiy[i], posiz[i],
	    velox[i], veloy[i], veloz[i]);
  }
  fprintf (fp, "%f %f %f\n", boxx, boxy, boxz);
}


