//My
#include "DetectorConstructionG2.hh"
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

G4VPhysicalVolume* DetectorConstructionG2::Construct()
{

  //ConstructField();
  
  G4double world_sizeX = 4.0*cm;
  G4double world_sizeY = world_sizeX;
  G4double world_sizeZ = 6.0*cm;

  //ARICH 
  G4double aerogel_body_Rin   = 0.0*mm;
  G4double aerogel_body_Rout  = 2/2.0*mm;
  G4double aerogel_body_X0 = 0.0*mm;
  G4double aerogel_body_Y0 = 0.0*mm;
  G4double aerogel_body_Z0 = 0.0*mm;

  //ARICH
  //G4double sensitive_Rin  = aerogel_body_Rin;
  //G4double sensitive_Rout = aerogel_body_Rout;
  //Tubs
  //G4double sensitive_Rin  = 0.0;
  //G4double sensitive_Rout = world_sizeX/2.0 - 1.0*mm;
  //Sphere
  G4double sensitive_Rout = world_sizeX/2.0 - 1.0*mm;
  G4double sensitive_Rin  = sensitive_Rout - 1.0*mm;

  G4double sensitive_X0 = 0.0;
  G4double sensitive_Y0 = 0.0;
  //G4double sensitive_Z0 = aerogel_body_Z0 + aerogel_body_sizeZ/2.0 + 0.1*cm;
  G4double sensitive_Z0 = 0.0;


  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;

  // 
  // Define World Volume
  //
  G4VSolid *world_solid = new G4Box("World",world_sizeX/2.0,world_sizeY/2.0,world_sizeZ/2.0);
  G4LogicalVolume *world_logical = new G4LogicalVolume(world_solid,Air,"World");
  G4VPhysicalVolume *world_physical = new G4PVPlacement(0,G4ThreeVector(),world_logical,"World",0,false,0);

  //
  // Aerogel
  //
  //G4VSolid *aerogel_body_solid = new G4Box("aerogel_body",aerogel_body_sizeX/2.0,aerogel_body_sizeY/2.0,aerogel_body_sizeZ/2.0);
  //G4VSolid *aerogel_body_solid = new G4Tubs("aerogel_body",aerogel_body_Rin,aerogel_body_Rout,aerogel_body_sizeZ/2.0,0, 360.0*deg);
  G4VSolid *aerogel_body_solid = new G4Sphere("aerogel_body",
					      aerogel_body_Rin,
					      aerogel_body_Rout,
					      0,
					      360.0*deg,
					      0,
					      180.0*deg);
  G4LogicalVolume *aerogel_body_logical = new G4LogicalVolume(aerogel_body_solid,AerogelA,"aerogel_body");
  Ta.setX(aerogel_body_X0);
  Ta.setY(aerogel_body_Y0);
  Ta.setZ(aerogel_body_Z0);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *aerogel_body_physical = new G4PVPlacement(Tr,                   //Transformation
							       aerogel_body_logical, //its logical volume				 
							       "aerogel_body",       //its name
							       world_logical,   //its mother  volume
							       false,                //no boolean operation
							       0);	             //copy number
  aerogel_body_physical->GetName();

  //
  // Sensitive volume
  //
  //G4VSolid *sensitive_solid = new G4Box("Sensitive", sensitive_sizeX/2.0, sensitive_sizeY/2.0, sensitive_sizeZ/2.0);
  //G4VSolid *sensitive_solid = new G4Tubs("aerogel_body",sensitive_Rin,sensitive_Rout,sensitive_sizeZ/2.0,0, 360.0*deg);
  G4VSolid *sensitive_solid = new G4Sphere("Sensitive",
					   sensitive_Rin,
					   sensitive_Rout,
					   0,
					   360.0*deg,
					   0,
					   180.0*deg);
  G4LogicalVolume *sensitive_logical = new G4LogicalVolume(sensitive_solid, Aluminum,"Sensitive");
  Ta.setX(sensitive_X0);
  Ta.setY(sensitive_Y0);
  Ta.setZ(sensitive_Z0);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *sensitive_physical = new G4PVPlacement(Tr,                //Transformation
							    sensitive_logical, //its logical volume				 
							    "Sensitive",       //its name
							    world_logical,     //its mother  volume
							    false,	       //no boolean operation
							    0);	               //copy number
  sensitive_physical->GetName();

  //
  // Set Visualization Attributes
  //
  //G4Color blue        = G4Color(0., 0., 1.);
  //G4Color green       = G4Color(0., 1., 0.);
  G4Color red         = G4Color(1., 0., 0.);
  G4Color white       = G4Color(1., 1., 1.);
  //G4Color cyan        = G4Color(0., 1., 1.);
  G4Color DircColor   = G4Color(0.0, 0.0, 1.0, 0.2);
  G4Color SensColor   = G4Color(0.0, 1.0, 1.0, 0.1);

  worldVisAtt->SetColor(white);
  worldVisAtt->SetVisibility(true);
  quartzVisAtt->SetColor(DircColor);
  quartzVisAtt->SetVisibility(true);
  
  sensitiveVisAtt->SetColor(red);
  sensitiveVisAtt->SetVisibility(true);
  absVisAtt->SetColor(SensColor);
  absVisAtt->SetVisibility(true);

  sensitive_logical->SetVisAttributes(sensitiveVisAtt);
  aerogel_body_logical->SetVisAttributes(quartzVisAtt);

  //world.logical->SetVisAttributes(worldVisAtt);
  
  //
  // Define Optical Borders
  //

  // Surface for killing photons at borders
  const G4int num1 = 2;
  G4double Ephoton[num1] = {1.5*eV, 5.8*eV};

  G4OpticalSurface* OpVolumeKillSurface =
    new G4OpticalSurface("VolumeKillSurface");
  OpVolumeKillSurface->SetType(dielectric_metal);
  OpVolumeKillSurface->SetFinish(polished);
  OpVolumeKillSurface->SetModel(glisur);


  G4double ReflectivityKill[num1] = {0., 0.};
  G4double EfficiencyKill[num1] = {1., 1.};
  G4MaterialPropertiesTable* VolumeKill = new G4MaterialPropertiesTable();
  VolumeKill->AddProperty("REFLECTIVITY", Ephoton, ReflectivityKill, num1);
  VolumeKill->AddProperty("EFFICIENCY",   Ephoton, EfficiencyKill,   num1);
  OpVolumeKillSurface->SetMaterialPropertiesTable(VolumeKill);
  new G4LogicalSkinSurface("SensitiveSurface", 
  			   sensitive_logical, OpVolumeKillSurface);
  


  // 
  // Sensitive detector definition
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  SensitiveDetector* aSD = new SensitiveDetector("fTOF");
  SDman->AddNewDetector(aSD);
  sensitive_logical->SetSensitiveDetector(aSD);
  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  G4double maxStep   = 0.1*mm;
  G4double maxLength = 2.0*m;
  G4double maxTime   = 20.0*ns; 
  G4double minEkin   = 1.0/100*MeV;
  G4double mionRang  = 0.01*mm;
  stepLimit = new G4UserLimits(maxStep,maxLength,maxTime,minEkin,mionRang);
  /*
  secA.logical->SetUserLimits(stepLimit);
  secB.logical->SetUserLimits(stepLimit);
  secC.logical->SetUserLimits(stepLimit);
  secWin.logical->SetUserLimits(stepLimit);
  */

  //G4GDMLParser parser;
  //parser.Write("CpFM.gdml", world.physical);

  return world_physical;
}
