//My
#include "DetectorConstructionG4.hh"
#include "SensitiveDetector.hh"
#include "MagneticField.hh"

//G4
#include "G4SystemOfUnits.hh"
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

using namespace CLHEP;

G4VPhysicalVolume* DetectorConstructionG4::Construct()
{

  ConstructField();
  
  G4double world_sizeX = 30.0*cm;
  G4double world_sizeY = 30.0*cm;
  G4double world_sizeZ = 50.0*cm;

  // G4Transform3D
  G4double volume1_sizeX  = 10.0*mm;
  G4double volume1_sizeY  =  5.0*mm;
  G4double volume1_sizeZ  =  1.0*mm;
  G4double volume1_X0     =  0.0*cm;
  G4double volume1_Y0     =  5.0*cm;
  G4double volume1_Z0     =  0.0*mm;
  G4double volume1_angle  = 45.0*deg;

  // G4RotationMatrix and G4ThreeVector
  G4double volume2_sizeX = 10.0*mm;
  G4double volume2_sizeY =  5.0*mm;
  G4double volume2_sizeZ =  1.0*mm;
  G4double volume2_X0    =  0.0*cm;
  G4double volume2_Y0    =  5.0*cm;
  G4double volume2_Z0    =  0.0*mm;
  G4double volume2_angle = 45.0*deg;

  // 
  // Define World Volume
  //
  G4VSolid *world_solid = new G4Box("World",world_sizeX/2.0,world_sizeY/2.0,world_sizeZ/2.0);
  G4LogicalVolume *world_logical = new G4LogicalVolume(world_solid,Air,"World");
  G4VPhysicalVolume *world_physical = new G4PVPlacement(0,G4ThreeVector(),world_logical,"World",0,false,0);

  //
  // volume1
  //
  G4VSolid *volume1_body_solid = new G4Box("volume1_body", volume1_sizeX/2.0, volume1_sizeY/2.0, volume1_sizeZ/2.0);
  G4LogicalVolume *volume1_body_logical = new G4LogicalVolume( volume1_body_solid, Air, "volume1_body_logical");
  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;
  Ta.setX(volume1_X0);
  Ta.setY(volume1_Y0);
  Ta.setZ(volume1_Z0);
  Ra.rotateZ(volume1_angle);
  Tr = G4Transform3D(Ra, Ta);
  //G4PVPlacement (const G4Transform3D &Transform3D, 
  //	 	   G4LogicalVolume *pCurrentLogical,
  //		   const G4String &pName,
  //		   G4LogicalVolume *pMotherLogical,
  //		   G4bool pMany, 
  //		   G4int pCopyNo,
  //		   G4bool pSurfChk=false)
  G4VPhysicalVolume *volume1_body_physical = new G4PVPlacement(Tr,                   //Transformation
							       volume1_body_logical, //its logical volume				 
							       "volume1_body",       //its name
							       world_logical,        //its mother  volume
							       false,                //no boolean operation
							       0);	             //copy number

  //
  // volume2
  //
  G4VSolid *volume2_body_solid = new G4Box("volume2_body", volume2_sizeX/2.0, volume2_sizeY/2.0, volume2_sizeZ/2.0);
  G4LogicalVolume *volume2_body_logical = new G4LogicalVolume( volume2_body_solid, Air, "volume2_body_logical");
  G4RotationMatrix *rotM = new G4RotationMatrix();
  G4ThreeVector trV;
  trV.setX(volume2_X0);
  trV.setY(volume2_Y0);
  trV.setZ(volume2_Z0);
  rotM->rotateZ(volume2_angle);
  //G4PVPlacement (G4RotationMatrix *pRot, 
  //		   const G4ThreeVector &tlate, 
  //		   G4LogicalVolume *pCurrentLogical, 
  //		   const G4String &pName, 
  //		   G4LogicalVolume *pMotherLogical, 
  //		   G4bool pMany, 
  //		   G4int pCopyNo, 
  //		   G4bool pSurfChk=false)
  G4VPhysicalVolume *volume2_body_physical = new G4PVPlacement(rotM,                 //Rotation matrix
							       trV,                  //Translation vector
							       volume2_body_logical, //its logical volume				 
							       "volume2_body",       //its name
							       world_logical,        //its mother  volume
							       false,                //no boolean operation
							       0);	             //copy number
  
  return world_physical;
}
