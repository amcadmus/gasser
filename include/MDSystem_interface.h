#ifndef __MDSystem_interface_h_wanghan__
#define __MDSystem_interface_h_wanghan__

#include "MDSystem.h"
#include "BondList_interface.h"
#include "AngleList_interface.h"
#include "BoxGeometry.h"
#include "NonBondedInteraction.h"
#include "MDTimer_interface.h"

#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"

using namespace RectangularBoxGeometry;

class MDSystem
{
  IndexType tmpNAtomType;
  // xdr variables
  XDRFILE *xdfile;
  matrix xdbox;
  rvec *xdx;
  float xdprec;
public:
  HostMDData   hdata;		/**< MD data on host */
  DeviceMDData ddata;		/**< MD data on device */
  DeviceMDData recoveredDdata;	/**< MD data on host. Used to hold the
				 * data recovered from the reshuffled
				 * system */
  RectangularBox box;		/**< box geometry */
  bool				hasBond;
  BondList			bdlist;		/**< bond list of the system */
  bool				hasAngle;
  AngleList			anglelist;	/**< angle list of the system */
  SystemNonBondedInteraction	nbInter;
public:
  MDSystem () ;
  ~MDSystem () ;
public:
  HostMDData & hostData () {return hdata;}
  const HostMDData & hostData () const {return hdata;}
  DeviceMDData & deviceData () {return ddata;}
  const DeviceMDData & deviceData () const {return ddata;}
public:
  /** 
   * Initialize system configuration from a .gro file.
   * 
   * @param configfile The name of the configuration file. It should be a .gro
   * file.
   * @param mapfile The atom name -> type mapping file.n
   * @param maxNumAtom Possible maximum number of atoms in the system. If value
   * 0 is given, then this value is supposed to be the same as the number of
   * atoms in the configuration file.
   */
  void initConfig (const char * configfile,
		   const char * mapfile,
		   const IndexType & maxNumAtom = 0);
  void setBoxSize (const ScalorType & x, const ScalorType & y, const ScalorType & z)
      {RectangularBoxGeometry::setBoxSize (x, y, z, &box);}

  void addNonBondedInteraction (const TypeType &i,
				const TypeType &j,
				const NonBondedInteractionParameter & param);
  void buildNonBondedInteraction ();
  ScalorType calMaxNBRcut ();
  void initBond ();
  void addBond (const IndexType & ii,
		const IndexType & jj,
		const BondInteractionParameter & param);
  void buildBond ();

  void initAngle ();
  void addAngle (const IndexType & ii,
		 const IndexType & jj,
		 const IndexType & kk,
		 const AngleInteractionParameter & param);
  void buildAngle ();
public:
  /** 
   * Update host data from device.
   * 
   * @param timer A timer monitors the time consumption.
   */
  void updateHost (MDTimer *timer=NULL);
  /** 
   * Initialize .xtc file output.
   * 
   * @param filename File name of the .xtc file.
   * @param prec Precision. A default value of 1000.f is assumed
   */
  void initWriteXtc (const char * filename, float prec = 1000.0f);
  /** 
   * Write host data to the .xtc file.
   * 
   * @param step Present step.
   * @param time Present time.
   * @param timer A timer monitors the time consumption.
   */
  void writeHostDataXtc (int step, float time, MDTimer *timer=NULL);
  /** 
   * Close the .xtc file.
   * 
   */
  void endWriteXtc ();
  
  void writeHostDataGro (const char * filename,
			 int step,
			 float time,
			 MDTimer *timer=NULL);
}
    ;





#endif
