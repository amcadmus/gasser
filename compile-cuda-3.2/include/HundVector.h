#ifndef __HundVector_h_wanghan__
#define __HundVector_h_wanghan__

class FailAllocateHundVector 
{
};

template <typename TYPE>
class HundVector
{
  TYPE * buff;
  unsigned dataSize;
  unsigned memSize;
public:
  HundVector ();
  HundVector (const HundVector & v);
  HundVector (const unsigned & size);
  HundVector (const unsigned & size, const TYPE & value);
public:
  void resize (const unsigned & size);
  void resize (const unsigned & size, const TYPE & value);
  HundVector & operator = (const HundVector & v);
};


template <typename TYPE>
HundVector<TYPE>::
HundVector ()
{
  buff = NULL;
  dataSize = 0;
  memSize = 0;
}

template <typename TYPE>
HundVector<TYPE>::
HundVector (const HundVector & v)
    : dataSize(v.dataSize), memSize(v.memSize)
{
  buff = NULL;
  buff = (TYPE *) realloc (buff, sizeof(TYPE) * memSize);
  if (buff == NULL){
    throw FailAllocateHundVector();
  }
  for (unsigned i = 0; i < dataSize; ++i){
    buff[i] = v.buff[i];
  }
}

template <typename TYPE>
HundVector<TYPE>::
HundVector (const unsigned & size)
    : dataSize(size), memSize(size)
{
  buff = NULL;
  buff = (TYPE *) realloc (buff, sizeof(TYPE) * memSize);
  if (buff == NULL){
    throw FailAllocateHundVector();
  }
}

template <typename TYPE>
HundVector<TYPE>::
HundVector (const unsigned & size, const TYPE & value)
    : dataSize(size), memSize(size)
{
  buff = NULL;
  buff = (TYPE *) realloc (buff, sizeof(TYPE) * memSize);
  if (buff == NULL){
    throw FailAllocateHundVector();
  }
  for (unsigned i = 0; i < size; ++i){
    buff[i] = value;
  }
}



#endif
