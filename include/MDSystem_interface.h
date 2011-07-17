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

/// The MDSystem contains data of all atom in a system.

/**
 * The MDSystem contains data of all atoms in a system, namely,
 * coordinates, velocities, forces, types, charges and masses as well
 * as the the box geometry. The box is presently rectangular.
 *
 * There are two copies of all the data, one is on the host (CPU
 * memory) and the other is on the device (GPU). MDSystem offers
 * several functions to transfer data between the host copy and device
 * copy.
 *
 * MDSystem also offers some input/output functions that faciliate the
 * recording of the MD trajactory (in .xtc file format) and the
 * configuration (in .gro file format). Remember all data
 * input/outputs are dealing with the host copy. One should first
 * transfer the device data to host before before any input/output of
 * the device data.
 * 
 */


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
  HostMDData   hdata;		/**< MD data on the host (host MD data) */
  DeviceMDData ddata;		/**< MD data on the device (device MD data)*/
  RectangularBox box;		/**< Simulation box */
public:
  MDSystem () ;
  ~MDSystem () ;
public:
  /** 
   * Get host MD data.
   * 
   * 
   * @return host MD data.
   */
  HostMDData & hostData () {return hdata;}
  /** 
   * Get host MD data.
   * 
   * 
   * @return host MD data.
   */
  const HostMDData & hostData () const {return hdata;}
  // DeviceMDData & deviceData () {return ddata;}
  // const DeviceMDData & deviceData () const {return ddata;}
public:
  /** 
   * Reshuffle the MDSystem.
   * 
   * @param indexTable Index table to reshuffle.
   * @param numAtom Number of atom in system, which is the same as the
   * length of the index table
   * @param timer Timer measuring the performance.
   */
  virtual void reshuffle (const IndexType * indexTable,
			  const IndexType & numAtom,
			  MDTimer * timer = NULL);
  /** 
   * Recover device MD data from the reshuffled system.
   * 
   * @param timer Timer measuring the performance.
   */
  void recoverDeviceData (MDTimer * timer = NULL);
public:
  /** 
   * Move all atoms out of the simulation box back.
   * 
   * @param timer Timer measuring the performance.
   */
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
  /** 
   * Initialize the type, mass, charge of all atoms.
   * 
   * @param sysTop 
   */
  void initTopology (const Topology::System & sysTop);
  /** 
   * Initialize device MD data from the host MD data.
   * 
   */
  void initDeviceData ();
  /** 
   * Set the rectangular simulation box of the system.
   * 
   * @param x Box size on x.
   * @param y Box size on y.
   * @param z Box size on z.
   */
  void setBoxSize (const ScalorType & x,
		   const ScalorType & y,
		   const ScalorType & z)
      {RectangularBoxGeometry::setBoxSize (x, y, z, &box);}
public:
  /** 
   * Update host MD data from device MD data.
   * 
   * @param timer A timer monitors the time consumption.
   */
  void updateHostFromDevice    (MDTimer *timer=NULL);
  /** 
   * Update host MD data from the recovered device MD data. (The device MD data
   * should be first recovered from reshuffled system by calling recoverDeviceData.)
   * 
   * @param timer A timer monitors the time consumption.
   */
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

  /** 
   * Write host data to .gro file.
   * 
   * @param filename file name.
   * @param step Present step.
   * @param time Present time.
   * @param timer A timer monitors the time consumption.
   */
  void writeHostDataGro (const char * filename,
			 int step,
			 float time,
			 MDTimer *timer=NULL);
}
    ;





#endif
