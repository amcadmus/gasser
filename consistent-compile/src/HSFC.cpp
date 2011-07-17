#define CPP_FILE
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>

// #include "MDException.h"
#include "common.h"
#include "HSFC.h"

template <typename TYPE >
class MatrixData
{
public:
  unsigned size;
public:
  std::vector<TYPE > data;
public:
  MatrixData ();
  MatrixData (const unsigned & size)
      { resize (size); }
  void resize (const unsigned & size);
  TYPE & ele (const unsigned & i, const unsigned & j, const unsigned & k)
      { return data[size * size * i + size * j + k]; }
  const TYPE & ele (const unsigned & i, const unsigned & j, const unsigned & k) const
      { return data[size * size * i + size * j + k]; }
};


template <typename TYPE >
MatrixData<TYPE>::
MatrixData()
{
  size = 0;
}

template <typename TYPE >
void MatrixData<TYPE>::
resize (const unsigned & nsize)
{
  size = nsize;
  data.resize (size * size * size);
}

class PRule 
{
  std::vector<unsigned > data;
public:
  PRule ();
public:
  void posiChange (const unsigned & relativeOrder);
  unsigned mapToOldIndex (const unsigned & index) {return data[index];}
};

PRule::
PRule()
{
  data.resize (8, 0);
  for (unsigned i = 0; i < 8; ++i){
    data[i] = i;
  }
}

void swap (unsigned & i, unsigned &j)
{
  unsigned tmp (i);
  i = j;
  j = tmp;
}

void PRule::
posiChange (const unsigned & posi)
{
  unsigned tmp;
  switch (posi%8){
  case 0:
      swap (data[1], data[7]);
      swap (data[2], data[4]);
      break;
  case 1:
      swap (data[2], data[6]);
      swap (data[3], data[7]);
      break;
  case 2:
      swap (data[2], data[6]);
      swap (data[3], data[7]);
      break;
  case 3:
      swap (data[0], data[2]);
      swap (data[1], data[3]);
      swap (data[4], data[6]);
      swap (data[5], data[7]);
      break;
  case 4:
      swap (data[0], data[2]);
      swap (data[1], data[3]);
      swap (data[4], data[6]);
      swap (data[5], data[7]);      
      break;
  case 5:
      swap (data[0], data[4]);
      swap (data[1], data[5]);
      break;
  case 6:
      swap (data[0], data[4]);
      swap (data[1], data[5]);
      break;
  case 7:
      swap (data[0], data[6]);
      swap (data[3], data[5]);
      break;
  }
}  
  

struct Mapping 
{
  MatrixData <unsigned > indexData;
  MatrixData <PRule >    ruleData;
  unsigned level;
public:
  unsigned simpleUnchangedFwd (const unsigned & i, const unsigned & j, const unsigned & k);
  void simpleUnchangedBkwd (const unsigned & index,
			    unsigned & i, unsigned & j, unsigned & k);
public:
  Mapping ();
  void levelUp ();
  void print ();
  void buildBackwardMapping (std::vector<std::vector<unsigned > > &bkmap);
};

Mapping::
Mapping ()
    : indexData (2), 
      ruleData  (2), 
      level (1)
{
  
  for (unsigned i = 0; i < 2; i ++){
    for (unsigned j = 0; j < 2; j ++){
      for (unsigned k = 0; k < 2; k ++){
	indexData.ele(i, j, k) = simpleUnchangedFwd (i, j, k);
      }
    }
  }
  for (unsigned i = 0; i < 2; i ++){
    for (unsigned j = 0; j < 2; j ++){
      for (unsigned k = 0; k < 2; k ++){
	ruleData.ele(i, j, k).posiChange (indexData.ele(i, j, k));
      }
    }
  }
}

unsigned Mapping::
simpleUnchangedFwd (const unsigned & i, const unsigned & j, const unsigned & k)
{
  switch ( (i << 2) + (j << 1) + k) {
  case 0: return 0;
  case 1: return 7;
  case 2: return 1;
  case 3: return 6;
  case 4: return 3;
  case 5: return 4;
  case 6: return 2;
  case 7: return 5;
  }
  return -1;
}

void Mapping::
simpleUnchangedBkwd (const unsigned & index,
		     unsigned & i, unsigned & j, unsigned & k)
{
  switch (index){
  case 0:
      i = 0;
      j = 0;
      k = 0;
      return;
  case 1:
      i = 0;
      j = 1;
      k = 0;
      return;
  case 2:
      i = 1;
      j = 1;
      k = 0;
      return;
  case 3:
      i = 1;
      j = 0;
      k = 0;
      return;
  case 4:
      i = 1;
      j = 0;
      k = 1;
      return;
  case 5:
      i = 1; 
      j = 1;
      k = 1;
      return;
  case 6:
      i = 0;
      j = 1;
      k = 1;
      return;
  case 7:
      i = 0;
      j = 0;
      k = 1;
      return;
  }
}
  

void Mapping::
levelUp ()
{
  MatrixData <unsigned > bkIndexData (indexData);
  MatrixData <PRule    > bkRuleData  (ruleData);
  unsigned bkLevel = level;

  level ++;
  unsigned presentSize = (1 << level);
  indexData.resize (presentSize);
  ruleData.resize  (presentSize);

  for (unsigned oldi = 0; oldi < (presentSize >> 1); ++oldi){
    for (unsigned oldj = 0; oldj < (presentSize >> 1); ++oldj){
      for (unsigned oldk = 0; oldk < (presentSize >> 1); ++oldk){
	for (unsigned relativeOrder = 0; relativeOrder < 8; ++ relativeOrder){
	  unsigned unchanged = bkRuleData.ele(oldi, oldj, oldk).mapToOldIndex(relativeOrder);
	  unsigned i , j, k;
	  simpleUnchangedBkwd (unchanged, i, j, k);
	  unsigned posii = (oldi << 1) + i;
	  unsigned posij = (oldj << 1) + j;
	  unsigned posik = (oldk << 1) + k;
	  indexData.ele(posii, posij, posik) = 
	      (bkIndexData.ele(oldi, oldj, oldk) << 3) + relativeOrder;
	  PRule tmpRule (bkRuleData.ele(oldi, oldj, oldk));
	  tmpRule.posiChange (relativeOrder);
	  ruleData.ele(posii, posij, posik) = tmpRule;
	}
      }
    }
  }
}

void Mapping::
print ()
{
  for (unsigned i = 0; i < indexData.size; ++i){
    for (unsigned j = 0; j < indexData.size; ++j){
      for (unsigned k = 0; k < indexData.size; ++k){
	printf ("%06d  ", indexData.ele(i, j, k));
      }
      printf ("\n");
    }
    printf ("\n");    
  }
}



void Mapping::
buildBackwardMapping (std::vector<std::vector<unsigned > > &bkmap)
{
  bkmap.resize(indexData.size * indexData.size * indexData.size);
  unsigned p = 0;
  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;
  for (; p < bkmap.size(); ++p){
    if (p != 0){
      for (int di = -1; di <= 1; di ++){
	for (int dj = -1; dj <= 1; dj ++){
	  for (int dk = -1; dk <= 1; dk ++){
	    int tryi = int(i) + di;
	    int tryj = int(j) + dj;
	    int tryk = int(k) + dk;
	    if (tryi >= 0 && tryi < indexData.size &&
		tryj >= 0 && tryj < indexData.size &&
		tryk >= 0 && tryk < indexData.size){
	      if (indexData.ele(tryi, tryj, tryk) == p){
		i = tryi;
		j = tryj;
		k = tryk;
		goto jumpout;
	      }
	    }
	  }
	}
      }
    }
    jumpout:
    std::vector<unsigned > tmp(3);
    tmp[0] = i;
    tmp[1] = j;
    tmp[2] = k;
    bkmap[p] = tmp;
  }
}
	

static unsigned Nbits (const unsigned & n)
{
  unsigned i = 0;
  while ((n >> i ) != 0) i++;
  return i;
}


static void HSFC3D (const unsigned & nx, const unsigned & ny, const unsigned & nz,
		    std::vector<std::vector<unsigned > > &bkmap)
{
  Mapping mp;
  unsigned bx = Nbits (nx);
  unsigned by = Nbits (ny);
  unsigned bz = Nbits (nz);
  
  unsigned max = bx;
  if (by > max) max = by;
  if (bz > max) max = bz;

  for (unsigned i = 0; i < max-1; ++i){
    mp.levelUp();
  }
  
  mp.buildBackwardMapping (bkmap);

  unsigned size = bkmap.size();
  for (unsigned i = 0; i < size; ++i){
    std::vector<unsigned > tmp = bkmap[i];
    if (tmp[0] >= nx || tmp[1] >= ny || tmp[2] >= nz){
      for (unsigned j = i; j < size-1; ++j){
	bkmap[j] = bkmap[j+1];
      }
      i --;
      size --;
    }
  }
  bkmap.resize(size);
}




void setHSFCMap1dto3d (HSFCMap1dto3d map,
		       IndexType index,
		       IndexType mp1, IndexType mp2, IndexType mp3)
{
  map[index * 3] = mp1;
  map[index * 3 + 1] = mp2;
  map[index * 3 + 2] = mp3;
}

void setHSFCMap3dto1d (HSFCMap3dto1d map,
		       IndexType nx, IndexType ny, IndexType nz,
		       IndexType i,  IndexType j,  IndexType k,
		       IndexType mp)
{
  map[i * ny * nz + j * nz + k] = mp;
}

void initHostHSFCMap1dto3d (HSFCMap1dto3d *map, 
			    IndexType nx, IndexType ny, IndexType nz )
{
  *map = (HSFCMap1dto3d) malloc (sizeof(IndexType) *nx*ny*nz * 3);
  if (*map == NULL){
    throw MDExcptFailedMallocOnHost("initHostHSFCMap1dto3d",
				    "map", sizeof(IndexType) *nx*ny*nz * 3);
  }
  
  std::vector<std::vector<unsigned > > bkmap;
  HSFC3D (nx, ny, nz, bkmap);

  for (unsigned i = 0; i < bkmap.size(); i++){
    setHSFCMap1dto3d (*map, i, bkmap[i][0], bkmap[i][1], bkmap[i][2]);
  }
}

void initHostHSFCMap3dto1d (HSFCMap3dto1d *map,
			    HSFCMap1dto3d map1to3,
			    IndexType nx, IndexType ny, IndexType nz )
{
  *map = (HSFCMap3dto1d) malloc (sizeof(IndexType) *nx*ny*nz * 3);
  if (*map == NULL){
    throw MDExcptFailedMallocOnHost("initHostHSFCMap3dto1d",
				    "map", sizeof(IndexType) *nx*ny*nz * 3);
  }
  for (unsigned ii = 0; ii < nx*ny*nz; ii += 1){
    unsigned i, j, k;
    map1dto3d (map1to3, ii, &i, &j, &k);
    setHSFCMap3dto1d (*map, nx, ny, nz, i, j, k, ii);
  }
}


