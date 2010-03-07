#ifndef __Parallel_CellList_h__
#define __Parallel_CellList_h__

#include "common.h"

namespace Parallel{

  class HostCellList 
  {
    HostVectorType frameLow;
    HostVectorType frameUp;
    HostIntVectorType numCell;
    HostVectorType numCelli;
    ScalorType rlist;

    IndexType * atomIndex;
    IndexType * numAtom;
    IndexType stride;

    IndexType * numNeighborCell;
    IndexType * neighborCellIndex;
    IndexType maxNumNeighborCell;
public:
    void reinit (const MDSystem  & sys,
		 const ScalorType & rlist,
		 const IndexType & NumThreadCell,
		 const IndexType & devideLevel = 1,
		 const BoxDirection_t & bdir = 7);
    void buildSubList (const IndexType & xIdLo,
		       const IndexType & xIdUp,
		       const IndexType & yIdLo,
		       const IndexType & yIdUp,
		       const IndexType & zIdLo,
		       const IndexType & zIdUp,
		       SubHostCellList & subList);
    void buildSubListLevel1 (const IndexType & xIdLo,
			     const IndexType & xIdUp,
			     const IndexType & yIdLo,
			     const IndexType & yIdUp,
			     const IndexType & zIdLo,
			     const IndexType & zIdUp,
			     SubHostCellList & subList);
  };

  class SubHostCellList 
  {
    HostCellList * motherList;
    IndexType * cellIndex;
    IndexType numCell;
    IndexType memSize;
    // IndexType * begin;
    // IndexType * end;
    class Iterator
    {
      IndexType * it;
  protected:
      Iterator (const IndexType * ptr) :it(ptr) {}
  public:
      Iterator ();
      Iterator (const Iterator & a) : it(a.it) {}
      Iterator & operator ++ () { ++it; return *this;}
      Iterator   operator ++ () const {Iterator tmpit(*this); it++; return tmpit;}
      Iterator & operator -- () { --it; return *this;}
      Iterator   operator -- () const {Iterator tmpit(*this); it--; return tmpit;}
      Iterator & operator =  (const Iterator & a) {return it =  a.it;}
      bool       operator == (const Iterator & a) {return it == a.it;}
      const IndexType & numCell () const;
      IndexType * atomIndexInCell () ;
      IndexType & numAtomInCell ();
    };
    friend class Iterator;

    Iterator begin () const {return begin;}
    Iterator end   () const {return end;}
public:
    void add (const SubHostCellList & a);
    void sub (const SubHostCellList & a);
    void merge (const SubHostCellList & a,
		const IndexType & shift);
    
  };
  
}

    


#endif

