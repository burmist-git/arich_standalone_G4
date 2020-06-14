#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

//G4
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"

class MagneticField;
class G4UserLimits;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  DetectorConstruction();
  ~DetectorConstruction();
   
public:
  
  G4VPhysicalVolume* Construct();
  void ConstructDetector();
  void ConstructField();

  //private:
protected:
  void DefineMaterials();
  G4UserLimits* stepLimit;        // pointer to user step limits

  //private:
protected:

  MagneticField *magField;

  void xmlarichdataUtil();

  // Various visibility attributes
  G4VisAttributes* worldVisAtt;
  G4VisAttributes* quartzVisAtt;
  G4VisAttributes* aerogelVisAtt;
  G4VisAttributes* sensitiveVisAtt;
  G4VisAttributes* pmtboxVisAtt;
  G4VisAttributes* absVisAtt;

  G4Material *Air;
  G4Material *C4F10;
  G4Material *Aerogel;
  G4Material *AerogelA;
  G4Material *AerogelB;
  G4Material *Aluminum;
  G4Material *AluminumMirr;
  G4Material *SiO2;
  G4Material* SiO2_cladd;
  G4Material* SiO2_coat; 
 
  /*
    WorldStruct world;
    SecStruct secA;
    SecStruct secB;
    SecStruct secC;
    SecStruct secAclad;
    SecStruct secBclad;
    SecStruct secCclad;
    SecStruct secAcoat;
    SecStruct secBcoat;
    SecStruct secCcoat;
    SecStruct secWin;
    SenDetStruct sensitive;
    Abs1Struct abs1;
    FiberStruct fiberCorr;
    FiberStruct fiberClad;
    FiberStruct fiberCoat;
    FiberStruct fiberBuff;
  */

  //LB need to be done stability tests
  //G4UserLimits* stepLimit;  // pointer to user step limits

};

#endif
