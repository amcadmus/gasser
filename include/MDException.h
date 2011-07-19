#ifndef __MDException_h_wanghan__
#define __MDException_h_wanghan__

// #include <string>
#include "string.h"
// #include "common.h"
#include <stdio.h>

class MDException;

class MDExcptCuda;

class MDExcptFile;
class MDExcptCannotOpenFile;
class MDExcptWrongFileFormat;

class MDExcptMemory;
class MDExcptFailedMallocOnHost;
class MDExcptFailedReallocOnHost;
class MDExcptExceedConstantMemLimit;

class MDExcptUndefinedAtomType;
class MDExcptUndefinedNBForceType;
class MDExcpt0AtomType;
class MDExcptUnbuiltNonBondedInteraction;
class MDExcptWrongNumberAtomDataTopology;


class MDExcptTooSmallCell;

// char delimitor[3] = {':', ' ', '\0'};




const char delimitor[4] = ": ";

class MDException 
{
public:
  virtual const char* what() const throw()
      {
	return "An exception happened";
      }
};

class MDExcptCuda 
{
public:
  virtual const char* what() const throw()
      {
	return "cuda exception";
      }
};

class MDExcpt0AtomType : public MDException
{
public:
  virtual const char* what() const throw()
      {
	return "Do not know how many types of atoms in the system";
      }
};

class MDExcptUndefinedAtomType : public MDException
{
  char message[MaxExceptionMsgLength];
public:
  MDExcptUndefinedAtomType (const char * fn) 
      { strncpy (message, "the type of atom ", MaxExceptionMsgLength);
	strncat (message, fn, MaxExceptionMsgLength);
	strncat (message, " is undefined", MaxExceptionMsgLength);}
  MDExcptUndefinedAtomType ( const char * description, const char * fn) 
      { strncpy (message, description, MaxExceptionMsgLength);
	strncat (message, delimitor, MaxExceptionMsgLength);
	strncat (message, "the type of atom ", MaxExceptionMsgLength);
	strncat (message, fn, MaxExceptionMsgLength); 
	strncat (message, " is undefined", MaxExceptionMsgLength);}
  virtual const char* what() const throw()
      {
	return message;
      }
};


class MDExcptUndefinedNBForceType : public MDException
{
  char message[MaxExceptionMsgLength];
public:
  MDExcptUndefinedNBForceType ( const char * description) 
      { strncpy (message, description, MaxExceptionMsgLength);
	strncat (message, delimitor, MaxExceptionMsgLength);
	strncat (message, "undefined non-bonded force type", MaxExceptionMsgLength);}
  virtual const char* what() const throw()
      {
	return message;
      }
};

class MDExcptUnbuiltNonBondedInteraction : public MDException
{
  char message[MaxExceptionMsgLength];
public:
  MDExcptUnbuiltNonBondedInteraction ( const char * description) 
      { strncpy (message, description, MaxExceptionMsgLength);
	strncat (message, delimitor, MaxExceptionMsgLength);
	strncat (message, "uninit non-boned interaction", MaxExceptionMsgLength);}
  virtual const char* what() const throw()
      {
	return message;
      }
};


class MDExcptFile :public MDException 
{
};

class MDExcptCannotOpenFile : public MDExcptFile
{
  char message[MaxExceptionMsgLength];
public:
  MDExcptCannotOpenFile (const char * fn) 
      { strncpy (message, "cannot open file ", MaxExceptionMsgLength);
	strncat (message, fn, MaxExceptionMsgLength); }
  MDExcptCannotOpenFile ( const char * description, const char * fn) 
      { strncpy (message, description, MaxExceptionMsgLength);
	strncat (message, delimitor, MaxExceptionMsgLength);
	strncat (message, "cannot open file ", MaxExceptionMsgLength);
	strncat (message, fn, MaxExceptionMsgLength); }
  virtual const char* what() const throw()
      {
	return message;
      }
};

class MDExcptWrongFileFormat : public MDExcptFile
{
  char message[MaxExceptionMsgLength];
public:
  MDExcptWrongFileFormat (const char * fn) 
      { strncpy (message, "wrong file format ", MaxExceptionMsgLength);
	strncat (message, fn, MaxExceptionMsgLength); }
  MDExcptWrongFileFormat ( const char * description, const char * fn) 
      { strncpy (message, description, MaxExceptionMsgLength);
	strncat (message, delimitor, MaxExceptionMsgLength);
	strncat (message, "wrong file format ", MaxExceptionMsgLength);
	strncat (message, fn, MaxExceptionMsgLength); }
  virtual const char* what() const throw()
      {
	return message;
      }
};

class MDExcptMemory : public MDException
{
public:
  virtual const char* what() const throw()
      {
	return "An unknow memory failure";
      }
};

class MDExcptFailedMallocOnHost : public MDExcptMemory
{
  char message[MaxExceptionMsgLength];
  size_t size;
public:
  MDExcptFailedMallocOnHost () {message[0] = '\0'; size = 0;}
  MDExcptFailedMallocOnHost (const char * name,
			     const size_t size) 
      { strncpy (message, "fail to malloc for ", MaxExceptionMsgLength);
	strncat (message, name, MaxExceptionMsgLength);
	sprintf (message, "%s with size %ld", message, size); }
  MDExcptFailedMallocOnHost (const char * description,
			     const char * name,
			     const size_t size) 
      { strncpy (message, description, MaxExceptionMsgLength);
	strncat (message, delimitor, MaxExceptionMsgLength);
	strncat (message, "fail to malloc for ", MaxExceptionMsgLength);
	strncat (message, name, MaxExceptionMsgLength);
	sprintf (message, "%s with size %ld", message, size); }
  virtual const char* what() const throw()
      {
	return message;
      }
};
  
  
class MDExcptFailedReallocOnHost : public MDExcptMemory
{
  char message[MaxExceptionMsgLength];
  size_t size;
public:
  MDExcptFailedReallocOnHost () {message[0] = '\0'; size = 0;}
  MDExcptFailedReallocOnHost (const char * name, const size_t size) 
      { strncpy (message, "fail to realloc for ", MaxExceptionMsgLength);
	strncat (message, name, MaxExceptionMsgLength);
	sprintf (message, "%s with size %ld", message, size); }
  MDExcptFailedReallocOnHost (const char * description,
			      const char * name,
			      const size_t size) 
      { strncpy (message, description, MaxExceptionMsgLength);
	strncat (message, delimitor, MaxExceptionMsgLength);
	strncat (message, "fail to realloc for ", MaxExceptionMsgLength);
	strncat (message, name, MaxExceptionMsgLength);
	sprintf (message, "%s with size %ld", message, size); }
  virtual const char* what() const throw()
      {
	return message;
      }
};
  
  
class MDExcptExceedConstantMemLimit : public MDExcptMemory
{
  char message[MaxExceptionMsgLength];
  size_t size;
public:
  MDExcptExceedConstantMemLimit () {message[0] = '\0'; size = 0;}
  MDExcptExceedConstantMemLimit (const char * name, const size_t size) 
      { strncpy (message, "exceed constant mem limit ", MaxExceptionMsgLength);
	strncat (message, name, MaxExceptionMsgLength);
	sprintf (message, "%s with size %ld", message, size); }
  MDExcptExceedConstantMemLimit (const char * description,
				 const char * name,
				 const size_t size) 
      { strncpy (message, description, MaxExceptionMsgLength);
	strncat (message, delimitor, MaxExceptionMsgLength);
	strncat (message, "exceed constant mem limit ", MaxExceptionMsgLength);
	strncat (message, name, MaxExceptionMsgLength);
	sprintf (message, "%s with size %ld", message, size); }
  virtual const char* what() const throw()
      {
	return message;
      }
};
  

class MDExcptTooSmallCell : public MDException
{
  char message[MaxExceptionMsgLength];
public:
  MDExcptTooSmallCell ()
      {
	message[0] ='\0';
    	strncat (message, "too small cell size", MaxExceptionMsgLength);
      }
  MDExcptTooSmallCell (const char * description,
		       const char * name,
		       const size_t size) 
      {
	message[0] ='\0';
	strncat (message, "in ", MaxExceptionMsgLength);
	strncat (message, description, MaxExceptionMsgLength);
	strncat (message, "too small cell size", MaxExceptionMsgLength);
      }
  virtual const char* what() const throw()
      {
	return message;
      }
};
  


class MDExcptWrongNumberAtomDataTopology : public MDException
{
public:
  virtual const char * what () const throw ()
      {
	return "Number of Atoms in topology and in data are not consistent";
      }
};
  



#endif
