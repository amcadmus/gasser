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

template <typename TYPE >
class MatrixData
{
public:
  unsigned stride;
  unsigned length;
public:
  std::vector<TYPE > data;
public:
  MatrixData ();
  MatrixData (const unsigned & stride, const unsigned & length)
      { resize (stride, length); }
  void resize (const unsigned & stride, const unsigned & length);
  TYPE & ele (const unsigned & i, const unsigned & j)
      { return data[stride * i + j]; }
  const TYPE & ele (const unsigned & i, const unsigned & j) const
      { return data[stride * i + j]; }
};


template <typename TYPE >
class Matrix
{
  unsigned stride;
  unsigned starti;
  unsigned startj;
  unsigned sizei;
  unsigned sizej;
  std::vector<TYPE > & refData;
public:
  Matrix (MatrixData<TYPE> & matData,
	  const unsigned & starti, const unsigned & startj,
	  const unsigned & sizei, const unsigned & sizej);
  void reinit (const unsigned & starti, const unsigned & startj,
	       const unsigned & sizei, const unsigned & sizej);
  TYPE & ele (const unsigned & i, const unsigned & j);
  const TYPE & ele (const unsigned & i, const unsigned & j) const;
};


template <typename TYPE >
MatrixData<TYPE>::
MatrixData()
{
  stride = length = 0;
}

template <typename TYPE >
void MatrixData<TYPE>::
resize (const unsigned & nstride, const unsigned & nlength)
{
  data.resize (nstride * nlength);
  stride = nstride;
  length = nlength;
}

template <typename TYPE >
Matrix<TYPE>::
Matrix (MatrixData<TYPE> & matData, 
	const unsigned & starti_, const unsigned & startj_,
	const unsigned & sizei_, const unsigned & sizej_)
    : refData (matData.data) ,
      starti(starti_), startj(startj_),
      sizei(sizei_), sizej(sizej_)
{
}

template <typename TYPE >
void Matrix<TYPE>::
reinit (const unsigned & starti_, const unsigned & startj_,
	const unsigned & sizei_, const unsigned & sizej_)
{
  stride = refData.stride;
  starti = starti_;
  startj = startj_;
  sizei = sizei_;
  sizej = sizej_;
}


template <typename TYPE >
inline TYPE & Matrix<TYPE>::
ele (const unsigned & i, const unsigned & j)
{
  return refData[(i+starti) * stride + j + startj];
}

template <typename TYPE >
inline const TYPE & Matrix<TYPE>::
ele (const unsigned & i, const unsigned & j) const
{
  return refData[(i+starti) * stride + j + startj];
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
  data.resize (4, 0);
  data[1] = 1;
  data[2] = 2;
  data[3] = 3;
}

void PRule::
posiChange (const unsigned & posi)
{
  unsigned tmp;
  switch (posi%4){
  case 0:
      tmp = (data[1]);
      data[1] = data[3];
      data[3] = tmp;
      break;
  case 1:
      break;
  case 2:
      break;
  case 3:
      tmp = (data[0]);
      data[0] = data[2];
      data[2] = tmp;
      break;
  }
}  
  

struct Mapping 
{
  MatrixData <unsigned > indexData;
  MatrixData <PRule >    ruleData;
  unsigned level;
public:
  unsigned simpleUnchangedFwd (const unsigned & i, const unsigned & j);
  void simpleUnchangedBkwd (const unsigned & index,
			    unsigned & i, unsigned & j);
public:
  Mapping ();
  void levelUp ();
  void print ();
  void buildBackwardMapping (std::vector<std::vector<unsigned > > &bkmap);
};

Mapping::
Mapping ()
    : indexData (2, 2), 
      ruleData  (2, 2), 
      level (1)
{
  indexData.ele(0, 0) = 0;
  indexData.ele(0, 1) = 1;
  indexData.ele(1, 1) = 2;
  indexData.ele(1, 0) = 3;
  for (unsigned i = 0; i < 2; i ++){
    for (unsigned j = 0; j < 2; j ++){
      ruleData.ele(i, j).posiChange (indexData.ele(i, j));
    }
  }
}

unsigned Mapping::
simpleUnchangedFwd (const unsigned & i, const unsigned & j)
{
  switch ( (i << 1) + j) {
  case 0: return 0;
  case 1: return 1;
  case 2: return 3;
  case 3: return 2;
  }
  return -1;
}

void Mapping::
simpleUnchangedBkwd (const unsigned & index,
		     unsigned & i, unsigned & j)
{
  switch (index){
  case 0:
      i = 0;
      j = 0;
      return;
  case 1:
      i = 0;
      j = 1;
      return;
  case 2:
      i = 1;
      j = 1;
      return;
  case 3:
      i = 1;
      j = 0;
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
  indexData.resize (presentSize, presentSize);
  ruleData.resize  (presentSize, presentSize);

  for (unsigned oldi = 0; oldi < (presentSize >> 1); ++oldi){
    for (unsigned oldj = 0; oldj < (presentSize >> 1); ++oldj){
      for (unsigned relativeOrder = 0; relativeOrder < 4; ++ relativeOrder){
	unsigned unchanged = bkRuleData.ele(oldi, oldj).mapToOldIndex(relativeOrder);
	unsigned i , j;
	simpleUnchangedBkwd (unchanged, i, j);
	unsigned posii = (oldi << 1) + i;
	unsigned posij = (oldj << 1) + j;
	indexData.ele(posii, posij) = 
	    (bkIndexData.ele(oldi, oldj) << 2) + relativeOrder;
	PRule tmpRule (bkRuleData.ele(oldi, oldj));
	tmpRule.posiChange (relativeOrder);
	ruleData.ele(posii, posij) = tmpRule;
      }
    }
  }
}

void Mapping::
print ()
{
  for (unsigned i = 0; i < indexData.length; ++i){
    for (unsigned j = 0; j < indexData.stride; ++j){
      printf ("%06d  ", indexData.ele(i, j));
    }
    printf ("\n");
  }
}



void Mapping::
buildBackwardMapping (std::vector<std::vector<unsigned > > &bkmap)
{
  bkmap.resize(indexData.stride * indexData.length);
  unsigned p = 0;
  unsigned i = 0;
  unsigned j = 0;
  for (; p < bkmap.size(); ++p){
    if (p != 0){
      for (int di = -1; di <= 1; di ++){
	for (int dj = -1; dj <= 1; dj ++){
	  int tryi = int(i) + di;
	  int tryj = int(j) + dj;
	  if (tryi >= 0 && tryi < indexData.length &&
	      tryj >= 0 && tryj < indexData.stride){
	    if (indexData.ele(tryi, tryj) == p){
	      i = tryi;
	      j = tryj;
	      goto jumpout;
	    }
	  }
	}
      }
    }
    jumpout:
    std::vector<unsigned > tmp(2);
    tmp[0] = i;
    tmp[1] = j;
    bkmap[p] = tmp;
  }
}
	
      

int main(int argc, char * argv[])
{
  Mapping mp;
  mp.levelUp();
  mp.levelUp();
  mp.levelUp();
  mp.levelUp();

  std::vector<std::vector<unsigned > > bkmap;
  mp.buildBackwardMapping (bkmap);
  
  // mp.print();

  for (unsigned i = 0; i < bkmap.size(); ++i){
    std::cout << bkmap[i][0] + 0.5 << " " << bkmap[i][1] + 0.5 << std::endl;
  }
  
}
