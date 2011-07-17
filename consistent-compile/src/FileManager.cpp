#define CPP_FILE

#include "FileManager.h"
#include "common.h"
#include <map>
#include <iostream>
#include <string.h>
#include <vector>
#include <algorithm>

struct TypeUnit
{
  TypeType type;
  ScalorType mass;
  ScalorType charge;
public:
  TypeUnit ()
      : type(0), mass (0), charge (0) {};
  TypeUnit (const TypeType & type_,
	    const ScalorType & mass_,
	    const ScalorType & charge_)
      : type(type_), mass (mass_), charge (charge_) {};
  bool operator == (const TypeUnit & u1)  {return this->type == u1.type;}
};


IndexType readAtomNameMapFile (const char * filename,
			       const IndexType numAtom,
			       const char * atomName,
			       TypeType * type,
			       ScalorType * mass,
			       ScalorType * charge)
{
  FILE * fp = fopen (filename, "r");
  if (fp == NULL) throw MDExcptCannotOpenFile ("readAtomNameMapFile", filename);
  char tmpname[StringSize];
  ScalorType tmpmass, tmpcharge;
  TypeType tmptype;
  int tag;
  std::vector<std::string> names;
  std::vector<TypeUnit   > data;

  char line [1024];
  while (fgets(line , 1024, fp) != NULL) {
    IndexType position = 0;
    while (isspace(line[position])) position ++;
    if (line[position] == '#') continue;
    tag = sscanf (line, "%s%d%f%f", 
		  tmpname, &tmptype, &tmpmass, &tmpcharge);
    if (tag == 4){
      TypeUnit tmpunit = TypeUnit (tmptype, tmpmass, tmpcharge);
      names.push_back (std::string(tmpname));
      data .push_back (tmpunit);
    }
  }
  
  std::vector<TypeUnit  > uniqueData (data);
  std::vector<TypeUnit  >::iterator new_end ;
  new_end = std::unique (uniqueData.begin(), uniqueData.end());
  uniqueData.erase (new_end, uniqueData.end());

  for (IndexType i = 0; i < names.size(); ++i){
    std::cout << "# "
	      << names[i] << "\t"
	      << data[i].type << " " 
	      << data[i].mass << " "
	      << data[i].charge << std::endl;
  }
  
  for (IndexType i = 0; i < numAtom; ++i){
    IndexType target = 0;
    for (; target < names.size(); ++target){
      if (names[target] == std::string(&atomName[i*StringSize])){
	type[i] = data[target].type;
	mass[i] = data[target].mass;
	charge[i] = data[target].charge;
	break;
      }
    }
    if (target == names.size()){
      throw MDExcptUndefinedAtomType("readAtomNameMapFile", &atomName[i*StringSize]);
    }
  }

  return uniqueData.size();
}


// IndexType readAtomNameMapFile (const char * filename,
// 			       const IndexType numAtom,
// 			       const char * atomName,
// 			       TypeType * type,
// 			       ScalorType * mass,
// 			       ScalorType * charge)
// {
//   FILE * fp = fopen (filename, "r");
//   if (fp == NULL) throw MDExcptCannotOpenFile ("readAtomNameMapFile", filename);
//   char tmpname[StringSize];
//   ScalorType tmpmass, tmpcharge;
//   TypeType tmptype;
//   int tag;
//   std::vector<std::string> names;
//   std::vector<TypeUnit   > data;
  
//   while ((tag=fscanf (fp, "%s%d%f%f", 
// 		      tmpname, &tmptype, &tmpmass, &tmpcharge)) != EOF){
//     if (tag != 4){
//       throw MDExcptWrongFileFormat ("readAtomNameMapFile", filename);
//     }
//     TypeUnit tmpunit = TypeUnit (tmptype, tmpmass, tmpcharge);
//     names.push_back (std::string(tmpname));
//     data .push_back (tmpunit);
//   }

//   std::vector<TypeUnit  > uniqueData (data);
//   std::vector<TypeUnit  >::iterator new_end ;
//   new_end = std::unique (uniqueData.begin(), uniqueData.end());
//   uniqueData.erase (new_end, uniqueData.end());
  
//   for (IndexType i = 0; i < numAtom; ++i){
//     IndexType target = 0;
//     for (; target < names.size(); ++target){
//       if (names[target] == std::string(&atomName[i*StringSize])){
// 	type[i] = data[target].type;
// 	mass[i] = data[target].mass;
// 	charge[i] = data[target].charge;
// 	break;
//       }
//     }
//     if (target == names.size()){
//       throw MDExcptUndefinedAtomType("readAtomNameMapFile", &atomName[i*StringSize]);
//     }
//   }

//   return uniqueData.size();
// }

