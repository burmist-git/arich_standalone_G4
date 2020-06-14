//My
#include "DetectorConstructionG1.hh"
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

G4VPhysicalVolume* DetectorConstructionG1::Construct()
{

  ConstructField();
  
  G4double world_sizeX = 4.0*m;
  G4double world_sizeY = 4.0*m;
  G4double world_sizeZ = 6.0*m;

  G4double c4f10_body_sizeX = 150*cm;
  G4double c4f10_body_sizeY = 120*2*cm;
  G4double c4f10_body_sizeZ = 93*cm;

  //Aerogel 
  //G4double aerogel_body_Rin = 40*cm;
  //G4double aerogel_body_Rout = 150*cm;
  G4double aerogel_body_sizeZ = 20.0*mm;
  G4double aerogel_body_Rin  = 441.0*mm;
  G4double aerogel_body_Rout = 1117.0*mm;
  G4double aerogel_body_X0 = 0.0*mm;
  G4double aerogel_body_Y0 = 0.0*mm;
  G4double aerogel_body_Z0 = 1668.0*mm + 10*mm;
  G4double aerogel_body_X0_A = aerogel_body_X0;
  G4double aerogel_body_Y0_A = aerogel_body_Y0;
  G4double aerogel_body_Z0_A = aerogel_body_Z0 + aerogel_body_sizeZ/2.0;
  G4double aerogel_body_X0_B = aerogel_body_X0;
  G4double aerogel_body_Y0_B = aerogel_body_Y0;
  G4double aerogel_body_Z0_B = aerogel_body_Z0_A + 1.0*mm + aerogel_body_sizeZ/2.0 + aerogel_body_sizeZ/2.0;

  //Mirror
  G4double flatMirror_lx = 429.6*mm;
  G4double flatMirror_ly = 199.0*mm;
  G4double flatMirror_lz = 2*mm;
  G4double flatMirror_X0 = 0.0*mm;
  G4double flatMirror_Y0 = 1136.0*mm - 10.0*mm;
  G4double flatMirror_Z0 = aerogel_body_Z0 + flatMirror_ly/2.0;
  G4double flatMirror_angle = 90.0*deg;

  //Sensitive
  G4double sensitive_Rin  = aerogel_body_Rin;
  G4double sensitive_Rout = aerogel_body_Rout;
  G4double sensitive_sizeZ = 2*mm;
  G4double sensitive_X0 = 0.0;
  G4double sensitive_Y0 = 0.0;
  G4double sensitive_Z0 = aerogel_body_Z0 + 20.0*cm;

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
  // Aerogel : A, B.
  //
  //G4VSolid *aerogel_body_solid = new G4Box("aerogel_body",aerogel_body_sizeX/2.0,aerogel_body_sizeY/2.0,aerogel_body_sizeZ/2.0);
  G4VSolid *aerogel_body_solid = new G4Tubs("aerogel_body",aerogel_body_Rin,aerogel_body_Rout,aerogel_body_sizeZ/2.0,0, 360.0*deg);
  G4LogicalVolume *aerogel_body_logical_A = new G4LogicalVolume(aerogel_body_solid,AerogelA,"aerogel_body_A");
  G4LogicalVolume *aerogel_body_logical_B = new G4LogicalVolume(aerogel_body_solid,AerogelB,"aerogel_body_B");
  //
  // Aerogel : A
  //
  Ta.setX(aerogel_body_X0_A);
  Ta.setY(aerogel_body_Y0_A);
  Ta.setZ(aerogel_body_Z0_A);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *aerogel_body_physical_A = new G4PVPlacement(Tr,                     //Transformation
								 aerogel_body_logical_A, //its logical volume				 
								 "aerogel_body_A",       //its name
								 world_logical,          //its mother  volume
								 false,                  //no boolean operation
								 0);	                 //copy number
  aerogel_body_physical_A->GetName();
  //
  // Aerogel : B
  //
  Ta.setX(aerogel_body_X0_B);
  Ta.setY(aerogel_body_Y0_B);
  Ta.setZ(aerogel_body_Z0_B);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *aerogel_body_physical_B = new G4PVPlacement(Tr,                     //Transformation
								 aerogel_body_logical_B, //its logical volume				 
								 "aerogel_body_B",       //its name
								 world_logical,          //its mother  volume
								 false,                  //no boolean operation
								 0);	                 //copy number
  aerogel_body_physical_B->GetName();


  //
  // Sensitive volume
  //
  //G4VSolid *sensitive_solid = new G4Box("Sensitive", sensitive_sizeX/2.0, sensitive_sizeY/2.0, sensitive_sizeZ/2.0);
  G4VSolid *sensitive_solid = new G4Tubs("aerogel_body",sensitive_Rin,sensitive_Rout,sensitive_sizeZ/2.0,0, 360.0*deg);
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

  //
  // Mirrors
  //
  G4VSolid *flatMirror_solid = new G4Box("flatMirror_solid",flatMirror_lx/2.0,flatMirror_ly/2.0,flatMirror_lz/2.0);
  G4LogicalVolume *flatMirror_logical = new G4LogicalVolume(flatMirror_solid,AluminumMirr,"flatMirror");
  //G4VPhysicalVolume *flatMirror_physical = new G4PVPlacement(0,G4ThreeVector(),world_logical,"flatMirror",0,false,0);
  Ta.setX(flatMirror_X0);
  Ta.setY(flatMirror_Y0);
  Ta.setZ(flatMirror_Z0);
  Ra.rotateX(flatMirror_angle);
  //Ra.rotateX(-1.0*deg);
  //Ra.rotateY(1.0*deg);
  //Ra.rotateZ(1.0*deg);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *flatMirror_physical = new G4PVPlacement(Tr,                 //Transformation
							     flatMirror_logical, //its logical volume				 
							     "flatMirror_body",  //its name
							     world_logical,      //its mother  volume
							     false,              //no boolean operation
							     0);	         //copy number
  flatMirror_physical->GetName();


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



  flatMirror_logical->SetVisAttributes(quartzVisAtt);
  sensitive_logical->SetVisAttributes(sensitiveVisAtt);
  aerogel_body_logical_A->SetVisAttributes(quartzVisAtt);
  aerogel_body_logical_B->SetVisAttributes(quartzVisAtt);
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
  



  // Define mirror surface
  const G4int num2 = 36;
  G4double EfficiencyMirrors[num2];
  G4double WaveLength[num2];
  G4double PhotonEnergy[num2];
  G4double MirrorReflectivity[num2];
  for (G4int i=0; i<num2; i++) {
    WaveLength[i] = (300 + i*10)*nanometer;
    PhotonEnergy[num2 - (i+1)] = twopi*hbarc/WaveLength[i];
    EfficiencyMirrors[i] = 0.0;
    MirrorReflectivity[i]=0.85;
  }
  /*
  G4double MirrorReflectivity[num2]=
    {0.87,0.88,0.885,0.89,0.895,0.9,0.905,0.91,0.915,0.92,0.923,0.9245,
     0.926,0.928,0.93,0.935,0.936,0.937,0.938,0.94,0.94,0.939,0.9382,
     0.938,0.937,0.937,0.936,0.935,0.934,0.932,0.93,0.928,0.926,0.924,
     0.922,0.92};
  */
  G4MaterialPropertiesTable* MirrorMPT = new G4MaterialPropertiesTable();
  MirrorMPT->AddProperty("RELECTIVITY", PhotonEnergy, MirrorReflectivity, num2);
  MirrorMPT->AddProperty("EFFICIENCY" , PhotonEnergy, EfficiencyMirrors,  num2);

  G4OpticalSurface* OpMirrorSurface = new G4OpticalSurface("MirrorSurface");
  OpMirrorSurface->SetType(dielectric_metal);
  OpMirrorSurface->SetFinish(polished);
  OpMirrorSurface->SetModel(glisur);

  new G4LogicalSkinSurface("MirrorSurfT",
			   flatMirror_logical, OpMirrorSurface);

  OpMirrorSurface->SetMaterialPropertiesTable(MirrorMPT);
  ///////////////////////////////////////
  
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
