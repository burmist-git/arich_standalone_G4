#ifndef DetectorConstructionG1_H
#define DetectorConstructionG1_H 1

//My
#include "DetectorConstruction.hh"

class DetectorConstructionG1 : public DetectorConstruction
{
public:
  
  DetectorConstructionG1() : DetectorConstruction()
  {
  }
     
  G4VPhysicalVolume* Construct();

};

#endif
