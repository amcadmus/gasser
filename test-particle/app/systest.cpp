#define HOST_CODE
#define CPP_FILE

#include "Parallel_MDSystem.h"
#include "NonBondedInteraction.h"
#include "GPU_Environment.h"

int main(int argc, char * argv[])
{
  Parallel::Environment env (&argc, &argv);
  int div[3];
  div[2] = env.numProc();
  div[1] = div[0] = 1;
  // div[0] = 2;
  // div[1] = 2;
  // div[2] = 2;
  env.init (div);

  GPU::Environment genv;
  genv.setDeviceId (0);
  
  IndexType nstep = 20;
  char * filename;  
  if (argc != 4){
    printf ("Usage:\n%s conf.gro nstep device\n", argv[0]);
    return 1;
  }
  if (argc != 1){
    nstep = atoi(argv[2]);
    filename = argv[1];
  }

  Parallel::MDSystem sys (env);
  
  Topology::System sysTop;
  Topology::Molecule mol;
  mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  LennardJones6_12Parameter ljparam;
  ljparam.reinit (1.f, 1.f, 0.f, 3.2f);
  mol.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam));
  sysTop.addMolecules (mol, sys.numAtomInGroFile(filename));

  sys.init (filename, sysTop);

  char name[1024];
  sprintf (name, "id%d.coord", env.myRank());
  sys.writeLocalData_SimpleFile (name);
  
  return 0;
}
