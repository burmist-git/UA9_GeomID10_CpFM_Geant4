#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

//My
#include "VolumeStructures.hh"

//G4
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"

class MagneticField;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  DetectorConstruction();
  ~DetectorConstruction();
   
public:
  
  G4VPhysicalVolume* Construct();
  void ConstructDetector();

private:
  void DefineMaterials();
  G4UserLimits* stepLimit;        // pointer to user step limits

private:

  //MagneticField *magField;

  // Various visibility attributes
  G4VisAttributes* worldVisAtt;
  G4VisAttributes* quartzVisAtt;
  G4VisAttributes* sensitiveVisAtt;
  G4VisAttributes* pmtboxVisAtt;
  G4VisAttributes* absVisAtt;

  WorldStruct world;
  SecStruct sec1;
  SecStruct sec2;
  SecStruct sec3;
  SenDetStruct sensitive;
  Abs1Struct abs1;

  Abs1Struct abs21;
  Abs1Struct abs22;
  Abs1Struct abs23;
  Abs1Struct abs24;

  FiberStruct fiberCorr;
  FiberStruct fiberClad;
  FiberStruct fiberCoat;
  FiberStruct fiberBuff;

  //LB need to be done stability tests
  //G4UserLimits* stepLimit;  // pointer to user step limits
};

#endif
