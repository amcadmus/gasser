class MDSystem;

#ifndef __MDSystem_interface_h_wanghan__
#define __MDSystem_interface_h_wanghan__

#include "MDSystem.h"
#include "BoxGeometry.h"
#include "MDTimer_interface.h"
#include "Interaction.h"
#include "Topology.h"
#include "Reshufflable.h"

#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"

using namespace RectangularBoxGeometry;

class MDSystem : public Reshufflable
{
  IndexType tmpNAtomType;
  // xdr variables
  XDRFILE *xdfile;
  matrix xdbox;
  rvec *xdx;
  float xdprec;
  IndexType * backMapTable;
  IndexType * backMapTableBuff;
  DeviceMDData recoveredDdata;	/**< MD data on device. Used to hold the
				 * data recovered from the reshuffled
				 * system */
  DeviceMDData bkDdata;		/**< MD data on device for reshuffle backup*/
public:
  HostMDData   hdata;		/**< MD data on host */
  DeviceMDData ddata;		/**< MD data on device */
  RectangularBox box;		/**< box geometry */
public:
  MDSystem () ;
  ~MDSystem () ;
public:
  HostMDData & hostData () {return hdata;}
  const HostMDData & hostData () const {return hdata;}
  // DeviceMDData & deviceData () {return ddata;}
  // const DeviceMDData & deviceData () const {return ddata;}
public:
  virtual void reshuffle (const IndexType * indexTable,
			  const IndexType & numAtom,
			  MDTimer * timer = NULL);
  void recoverDeviceData (MDTimer * timer = NULL);
public:
  void normalizeDeviceData (MDTimer * timer = NULL);
public:
  /** 
   * Initialize system configuration from a .gro file.
   * 
   * @param configfile The name of the configuration file. It should be a .gro
   * file.
   * @param maxNumAtom Possible maximum number of atoms in the system. If value
   * 0 is given, then this value is supposed to be the same as the number of
   * atoms in the configuration file.
   */
  void initConfig (const char * configfile,
		   // const char * mapfile,
		   const IndexType & maxNumAtom = 0);
  void initTopology (const Topology::System & sysTop);
  void initDeviceData ();
  void setBoxSize (const ScalorType & x, const ScalorType & y, const ScalorType & z)
      {RectangularBoxGeometry::setBoxSize (x, y, z, &box);}
public:
  /** 
   * Update host data from device.
   * 
   * @param timer A timer monitors the time consumption.
   */
  void updateHostFromDevice    (MDTimer *timer=NULL);
  void updateHostFromRecovered (MDTimer *timer=NULL);
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
