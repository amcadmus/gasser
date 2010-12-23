#define CPP_FILE

#include "Topology.h"

class MDExcptWrongAtomIndex : public MDException
{
public:
  virtual const char * what () const throw ()
      {
	return "wrong Atom Index";
      }
};

class MDExcptDoubleSpecifyBondInteraction : public MDException
{
public:
  virtual const char * what () const throw ()
      {
	return "the bond parameter should only be specified once";
      }
};

class MDExcptDoubleSpecifyAngleInteraction : public MDException
{
public:
  virtual const char * what () const throw ()
      {
	return "the angle parameter should only be specified once";
      }
};


Topology::Atom::
Atom()
    : mass (1.), charge (0.), type(0)
{
  name[0] = '\0';
}

Topology::Atom::
Atom (const Atom & a)
    : mass(a.mass), charge(a.charge), type(a.type)
{
}

Topology::Atom::
Atom (const ScalorType & mass_,
      const ScalorType & charge_,
      const TypeType & type_)
    : mass (mass_), charge(charge_), type(type_)
{
}

void Topology::Atom::
setProperty  (const ScalorType & mass_,
	      const ScalorType & charge_,
	      const TypeType & type_)
{
  mass = mass_;
  charge = charge_;
  type = type_;
}

Topology::NonBondedInteraction::
NonBondedInteraction (const TypeType & atomType0_,
		      const TypeType & atomType1_,
		      const NonBondedInteractionParameter & p)
    : atomType0(atomType0_), atomType1(atomType1_), ptr_nbInter(&p)
{
}

void Topology::NonBondedInteraction::
specifyInteraction(const TypeType & atomType0_,
		   const TypeType & atomType1_,
		   const NonBondedInteractionParameter & p)
{
  atomType0 = atomType0_;
  atomType1 = atomType1_;
  ptr_nbInter = &p;
}

// Topology::NonBondedInteraction::
// NonBondedInteraction (const TypeType & atomType0_,
// 		      const TypeType & atomType1_,
// 		      const NonBondedInteractionParameter & p)
//     : atomType0(atomType0_), atomType1(atomType1_), type(p.type())
// {
//   for (unsigned i = 0; i < p.numParam(); ++i){
//     paramArray.push_back( p.c_ptr()[i] );
//   }
  
//   rcut = p.rcut();
// }

// void Topology::NonBondedInteraction::
// specifyInteraction(const TypeType & atomType0_,
// 		   const TypeType & atomType1_,
// 		   const NonBondedInteractionParameter & p)
// {
//   atomType0 = atomType0_;
//   atomType1 = atomType1_;
//   type = p.type();
//   paramArray.clear();
//   for (unsigned i = 0; i < p.numParam(); ++i){
//     paramArray.push_back( p.c_ptr()[i] );
//   }
//   rcut = p.rcut();
// }

Topology::Exclusion::
Exclusion (const IndexType & atom0_,
	   const IndexType & atom1_)
    : atom0 (atom0_), atom1(atom1_)
{
}

Topology::Bond::
Bond (const IndexType & atom0_,
      const IndexType & atom1_,
      const BondInteractionParameter & p)
    :atom0(atom0_), atom1(atom1_), type(p.type())
{
  for (unsigned i = 0; i < p.numParam(); ++i){
    paramArray.push_back( p.c_ptr()[i] );
  }
}

void Topology::Bond::
specifyInteraction (const IndexType & atom0_,
		    const IndexType & atom1_,
		    const BondInteractionParameter & p)
{
  atom0 = atom0_;
  atom1 = atom1_;
  type = p.type();
  paramArray.clear();
  for (unsigned i = 0; i < p.numParam(); ++i){
    paramArray.push_back( p.c_ptr()[i] );
  }
}


Topology::Angle::
Angle (const IndexType & atom0_,
       const IndexType & atom1_,
       const IndexType & atom2_,
       const AngleInteractionParameter & p)
    :edge0(atom0_), center(atom1_), edge1(atom2_), type(p.type())
{
  for (unsigned i = 0; i < p.numParam(); ++i){
    paramArray.push_back( p.c_ptr()[i] );
  }
}

void Topology::Angle::
specifyInteraction (const IndexType & atom0_,
		    const IndexType & atom1_,
		    const IndexType & atom2_,
		    const AngleInteractionParameter & p)
{
  edge0 = atom0_;
  center = atom1_;
  edge1 = atom2_;
  type = p.type();
  paramArray.clear();
  for (unsigned i = 0; i < p.numParam(); ++i){
    paramArray.push_back( p.c_ptr()[i] );
  }
}      


Topology::Molecule::
Molecule ()
{
  name[0] = '\0';
}

void Topology::Molecule::
pushAtom (const Atom & a)
{
  atoms.push_back(a);
}

void Topology::System::
addNonBondedInteraction (const NonBondedInteraction & nb)
{
  nonBondedInteractions.push_back(nb);
  // printf ("#%f\n", nb.energyCorrection(1.8));
}

void Topology::Molecule::
addBond (const Bond & bd)
{
  if (bd.atom0 >= atoms.size() ||
      bd.atom1 >= atoms.size()){
    throw MDExcptWrongAtomIndex ();
  }
  bonds.push_back(bd);
}

void Topology::Molecule::
addAngle (const Angle & ag)
{
  if (ag.edge0 >= atoms.size() ||
      ag.edge1 >= atoms.size() ||
      ag.center>= atoms.size()){
    throw MDExcptWrongAtomIndex ();
  }
  angles.push_back(ag);
}

void Topology::Molecule::
addExclusion (const Exclusion & ex)
{
  if (ex.atom0 >= atoms.size() ||
      ex.atom1 >= atoms.size()){
    throw MDExcptWrongAtomIndex ();
  }
  exclusions.push_back(ex);
}

void Topology::Molecule::
clear()
{
  atoms.clear();
  bonds.clear();
  angles.clear();
}


Topology::System::
System ()
{
  name[0] = '\0';
  numFreedom = 0;
}


void Topology::System::
addMolecules (const Molecule & mol,
	      const IndexType & number)
{
  molecules.push_back(mol);
  numbers.push_back(number);
  if (indexShift.size() == 0) indexShift.push_back(0);
  indexShift.push_back (indexShift.back() + number * mol.atoms.size());
  numFreedom += mol.size() * number * 3;
}

void Topology::System::
clear()
{
  name[0] = '\0';
  numFreedom = 0;
  molecules.clear();
  numbers.clear();
  indexShift.clear();
}

void Topology::System::
calMolTopPosition (const IndexType & globalIndex,
		   IndexType & molIndex,
		   IndexType & atomIndex) const
{
  if (globalIndex >= indexShift.back ()){
    throw MDExcptTopology ("wrong global index");
  }
  molIndex = 0;
  while (globalIndex >= indexShift[molIndex+1]) molIndex++;
  atomIndex = (globalIndex - indexShift[molIndex]) % (molecules[molIndex].size());
}

