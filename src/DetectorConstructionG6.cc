//My
#include "DetectorConstructionG6.hh"
#include "SensitiveDetector.hh"
#include "MagneticField.hh"
#include "libxmlarichdata.h"

//G4
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Paraboloid.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Color.hh"
#include "G4TwoVector.hh"
#include "G4SDManager.hh"
#include "globals.hh"
//magnetic field
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4UserLimits.hh"
//GDML
//#include <G4GDMLParser.hh>

//root 
#include "TMath.h"

#include "G4SystemOfUnits.hh"

G4VPhysicalVolume* DetectorConstructionG6::Construct()
{

  xmlarichdata::quartz_Babar_DIRC_polishing_quality_bounds();
  
  ConstructField();  

  G4double world_sizeX = 10.0*cm;
  G4double world_sizeY = 10.0*cm;
  G4double world_sizeZ = 10.0*cm;

  G4double feb_alcooling_box1_sizeX = 20.5*mm;
  G4double feb_alcooling_box1_sizeY = 20.5*mm;
  G4double feb_alcooling_box1_sizeZ = 1.0*mm;
  
  G4double feb_alcooling_box2_sizeX = 32.4*mm;
  G4double feb_alcooling_box2_sizeY = 32.4*mm;
  G4double feb_alcooling_box2_sizeZ = feb_alcooling_box1_sizeZ;

  G4double feb_alcooling_box3_sizeX = 14.2*mm;
  G4double feb_alcooling_box3_sizeY = 31.8*mm;
  G4double feb_alcooling_box3_sizeZ = 0.5*mm;

  G4double feb_alcooling_box1_X0 = feb_alcooling_box2_sizeX/2.0 + feb_alcooling_box1_sizeX/2.0;
  G4double feb_alcooling_box1_Y0 = feb_alcooling_box2_sizeY/2.0 + feb_alcooling_box1_sizeY/2.0;
  G4double feb_alcooling_box1_Z0 = 0.0*mm;

  G4double feb_alcooling_box2_X0 = 0.0*mm;
  G4double feb_alcooling_box2_Y0 = 0.0*mm;
  G4double feb_alcooling_box2_Z0 = 0.0*mm;
  
  G4double feb_alcooling_box3_X0 = 27.0/TMath::Sqrt(2.0)*mm;
  G4double feb_alcooling_box3_Y0 = 27.0/TMath::Sqrt(2.0)*mm;
  G4double feb_alcooling_box3_Z0 = feb_alcooling_box1_sizeZ/2.0 + feb_alcooling_box3_sizeZ/2.0;
  G4double feb_alcooling_box3_angle = 45.0*deg;

  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;

  // 
  // Define World Volume
  //
  G4VSolid *world_solid = new G4Box("World",world_sizeX/2.0,world_sizeY/2.0,world_sizeZ/2.0);
  G4LogicalVolume *world_logical = new G4LogicalVolume(world_solid,Air,"World");
  G4VPhysicalVolume *world_physical = new G4PVPlacement(0,G4ThreeVector(),world_logical,"World",0,false,0);

  G4VSolid *feb_alcooling_box1_solid = new G4Box("feb_alcooling_box1_solid", feb_alcooling_box1_sizeX/2.0, feb_alcooling_box1_sizeY/2.0, feb_alcooling_box1_sizeZ/2.0);
  G4VSolid *feb_alcooling_box2_solid = new G4Box("feb_alcooling_box2_solid", feb_alcooling_box2_sizeX/2.0, feb_alcooling_box2_sizeY/2.0, feb_alcooling_box2_sizeZ/2.0);
  G4VSolid *feb_alcooling_box3_solid = new G4Box("feb_alcooling_box3_solid", feb_alcooling_box3_sizeX/2.0, feb_alcooling_box3_sizeY/2.0, feb_alcooling_box3_sizeZ/2.0);

  G4LogicalVolume *feb_alcooling_box1_logical = new G4LogicalVolume( feb_alcooling_box1_solid, Aluminum, "feb_alcooling_box1");
  G4LogicalVolume *feb_alcooling_box2_logical = new G4LogicalVolume( feb_alcooling_box2_solid, Aluminum, "feb_alcooling_box2");
  G4LogicalVolume *feb_alcooling_box3_logical = new G4LogicalVolume( feb_alcooling_box3_solid, Aluminum, "feb_alcooling_box3");

  //
  //Box 2
  //
  Ta.setX(feb_alcooling_box2_X0);
  Ta.setY(feb_alcooling_box2_Y0);
  Ta.setZ(feb_alcooling_box2_Z0);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *feb_alcooling_box2_physical = new G4PVPlacement(Tr,                         //Transformation
								     feb_alcooling_box2_logical, //its logical volume                            
								     "feb_alcooling_box2",       //its name
								     world_logical,              //its mother  volume
								     false,                      //no boolean operation
								     0);                         //copy number
  //
  //Box 1 A
  //
  Ta.setX(feb_alcooling_box1_X0);
  Ta.setY(feb_alcooling_box1_Y0);
  Ta.setZ(feb_alcooling_box1_Z0);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *feb_alcooling_box1_physical_A = new G4PVPlacement(Tr,                         //Transformation
								       feb_alcooling_box1_logical, //its logical volume                            
								       "feb_alcooling_box1_A",     //its name
								       world_logical,              //its mother  volume
								       false,                      //no boolean operation
								       0);                         //copy number

  //
  //Box 1 B
  //
  Ta.setX(-feb_alcooling_box1_X0);
  Ta.setY(-feb_alcooling_box1_Y0);
  Ta.setZ( feb_alcooling_box1_Z0);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *feb_alcooling_box1_physical_B = new G4PVPlacement(Tr,                         //Transformation
								       feb_alcooling_box1_logical, //its logical volume                            
								       "feb_alcooling_box1_B",     //its name
								       world_logical,              //its mother  volume
								       false,                      //no boolean operation
								       0);                         //copy number

  //
  //Box 3 A
  //
  Ta.setX(feb_alcooling_box3_X0);
  Ta.setY(feb_alcooling_box3_Y0);
  Ta.setZ(feb_alcooling_box3_Z0);
  Ra.rotateZ(-feb_alcooling_box3_angle);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *feb_alcooling_box3_physical_A = new G4PVPlacement(Tr,                         //Transformation
								       feb_alcooling_box3_logical, //its logical volume                            
								       "feb_alcooling_box3_A",     //its name
								       world_logical,              //its mother  volume
								       false,                      //no boolean operation
								       0);                         //copy number

  Ra.rotateZ(feb_alcooling_box3_angle);
  //
  //Box 3 B
  //
  Ta.setX(-feb_alcooling_box3_X0);
  Ta.setY(-feb_alcooling_box3_Y0);
  Ta.setZ(feb_alcooling_box3_Z0);
  Ra.rotateZ(-feb_alcooling_box3_angle);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *feb_alcooling_box3_physical_B = new G4PVPlacement(Tr,                         //Transformation
								       feb_alcooling_box3_logical, //its logical volume                            
								       "feb_alcooling_box3_B",     //its name
								       world_logical,              //its mother  volume
								       false,                      //no boolean operation
								       0);                         //copy number
  Ra.rotateZ(feb_alcooling_box3_angle);
  
  G4double maxStep   = 0.1*mm;
  G4double maxLength = 2.0*m;
  G4double maxTime   = 20.0*ns; 
  G4double minEkin   = 1.0/100*MeV;
  G4double mionRang  = 0.01*mm;
  stepLimit = new G4UserLimits(maxStep,maxLength,maxTime,minEkin,mionRang);

  return world_physical;
}
