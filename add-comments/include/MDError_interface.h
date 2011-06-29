#ifndef __DMError_interface_h_wanghan__
#define __DMError_interface_h_wanghan__

#include "common.h"
#include "error.h"

class MDError
{
  mdError_t he;
  char * getErrorString (mdError_t err);
  IndexType *hindex;
  ScalorType *hscalor;
public:
  IndexType  * ptr_dindex;
  ScalorType * ptr_dscalor;
  IndexType getIndex ();
  ScalorType getScalor ();
public:
  mdError_t * ptr_de;
public:
  MDError();
  ~MDError();
  void updateHost ();
  void check (const char * msg);
}
    ;


#endif
