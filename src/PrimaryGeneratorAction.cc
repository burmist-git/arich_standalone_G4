//my
#include "PrimaryGeneratorAction.hh"
#include "libxmlarichdata.h"

//G4
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"

#include "G4SystemOfUnits.hh"

//root
#include "TMath.h"
#include "TVector2.h"

using namespace CLHEP;

PrimaryGeneratorAction::PrimaryGeneratorAction() :
  _particleGun(0),
  _particleName("pi+"),
  _particleMomentum(3.0*GeV),
  _PhiAngle(0.0*deg),
  _ThetaAngle(0.0*deg),
  _X0(0.0*cm),
  _Y0(0.0*cm),
  _Z0(0.0*cm),
  _singlePhoton(false)
{
  _particleGun = new G4ParticleGun(1);  
  _BunchXID = 0;

  //backGen = new backgroundGen();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete _particleGun;
  //delete backGen;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  //xmlarichdata::hapdQE();

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle;
  // Correct for center of bar
  G4double xInit, yInit, zInit;
  G4double dX, dY, dZ;
  G4double Ekin, m;
  //G4int pdgID;
  //G4int i;
  //i = 0;

  //_particleMomentum = 7000.0*GeV;
  //_particleMomentum = 20.0*GeV;

  //xInit = 0.0*cm;
  //yInit = 0.0*cm;
  //zInit = -2.0*cm;
  xInit = _X0*cm;
  yInit = _Y0*cm;
  zInit = _Z0*cm;

  _BunchXID++;
  particle = particleTable->FindParticle(_particleName);
  m = particle->GetPDGMass();
  Ekin = (TMath::Sqrt(_particleMomentum*_particleMomentum + m*m) - m);

  //G4cout<<"_particleMomentum = "<<_particleMomentum<<G4endl
  //	<<"_ThetaAngle       = "<<_ThetaAngle*180.0/TMath::Pi()<<G4endl
  //	<<"_PhiAngle         = "<<_PhiAngle*180.0/TMath::Pi()<<G4endl;

  //_ThetaAngle = _ThetaAngle + (-1 + 2*G4UniformRand())*2*TMath::Pi()/360;//mearing of one degree
  //_PhiAngle   = _PhiAngle   + (-1 + 2*G4UniformRand())*2*TMath::Pi()/360;//mearing of one degree

  dX =  TMath::Sin(_ThetaAngle)*TMath::Cos(_PhiAngle);
  dY =  TMath::Sin(_ThetaAngle)*TMath::Sin(_PhiAngle);
  dZ =  TMath::Cos(_ThetaAngle);

  G4ThreeVector dir(dX, dY, dZ);
  _particleGun->SetParticleDefinition(particle);
  _particleGun->SetParticleMomentumDirection(dir);
  if(_particleName == "opticalphoton"){
    ////////////////////////
    //Generate single photon
    //Please note in new geant polarisation of the photon does not need to be defined
    //https://svn.lal.in2p3.fr/projects/GRED/Geant4Sim/fTOF_nPhot/fTOF_V0.0.4_svn_PvsTheta/src/fTOF_PrimaryGeneratorAction.cc
    //_particleName = "opticalphoton";
    //G4double wavelength = 410*nm;
    //_particleMomentum = twopi*hbarc/wavelength;
    ///////////////////////
    G4double wavelength = _particleMomentum*nm;
    //G4double wavelength = 420*nm;
    _particleGun->SetParticleEnergy(twopi*hbarc/wavelength);
    //Need to assign a polarization vector to optical photons
    //Polarisation of the Cherenkov photon is always perpendicular to cherenkov photon
    //and contained by the plane defined by the photon and particle 
    //So both particle and photon define polarisation vector
    //https://en.wikipedia.org/wiki/Cross_product --> see : The vector triple product
    //Explenatory picture : leonid@nb-gred02:/home/leonid/Dropbox/Thesis/tesina-icp.ps.gz page (9, 78) 
    //Formula :
    //pv -> polarisation vector of the photon (please note with this formula the vector is not unit vector)
    //gv -> unit vector collinear with photon direction 
    //mv -> unit vector collinear with particle direction 
    //Formula : pv = [gv x [ gv x mv]] = gv*(gv*mv) - mv
    //For single photon generator we use photon along z axis and the particle direction is irrelevant. 
    //So we arbitrary choose polarisation vector in y-direction.
    //Please _NOTE_ : in case direction of generated photon different the polarisation have to be defined in proper way.
    //(see the Formula : pv = [gv x [ gv x mv]] = gv*(gv*mv) - mv)
    G4double phi = (G4UniformRand() - 0.5)*twopi;
    TVector2 v2;
    v2.SetMagPhi(1.0,phi);
    G4double pol_x = v2.X();
    G4double pol_y = v2.Y();
    G4double pol_z = 0.0;
    G4ThreeVector polarization(pol_x, pol_y, pol_z);
    _particleGun->SetParticlePolarization(polarization);
    //G4cout<<"hbarc "<<hbarc<<G4endl;
  }
  else{
    _particleGun->SetParticleEnergy(Ekin);
  }
  _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
  //for(i = 0;i<10;i++)
  _particleGun->GeneratePrimaryVertex(anEvent);

}

G4int PrimaryGeneratorAction::GenFlatInt(G4int iMin,G4int iMax){
  G4int val;
  val = (G4int)floor((iMax - iMin + 1)*G4UniformRand() + iMin);
  return val;
}

void PrimaryGeneratorAction::generateThetaAndPhi(){
  _PhiAngle = G4UniformRand()*2*TMath::Pi();
  _ThetaAngle = TMath::Pi() - genCos2dist();
}

G4double PrimaryGeneratorAction::genCos2dist(){
  G4double theta = -999.0;//deg 
  G4double x = -999.0;
  G4double y = -999.0;
  while(theta==-999.0){
    x = G4UniformRand()*(70.0*TMath::Pi()/180.0); //rad
    y = G4UniformRand();
    if(TMath::Power(TMath::Cos(x),1.85)>y){
      theta = x;
    }
  }  
  return theta;
}
