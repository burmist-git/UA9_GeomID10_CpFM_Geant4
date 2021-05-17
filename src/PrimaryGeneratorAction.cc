//my
#include "PrimaryGeneratorAction.hh"
#include "VolumeStructures.hh"

//G4
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"

//root
#include "TMath.h"

PrimaryGeneratorAction::PrimaryGeneratorAction() :
  _particleGun(0),
  _particleName("pi+"),
  _particleMomentum(3.0*GeV),
  _PhiAngle(0.0*deg),
  _ThetaAngle(0.0*deg),
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

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle;
  // Correct for center of bar
  G4double xInit, yInit, zInit;
  G4double dX, dY, dZ;
  G4double Ekin, m;
  G4int pdgID;
  G4int i;

  //xInit = (G4UniformRand() - 0.5)*20.0*cm; 

  _particleMomentum = 1.0*MeV;
  //_particleMomentum = 7000.0*GeV;
  //_particleMomentum = 446.0*MeV;

  //with 90 deg anle - simple box
  //xInit =-4.0*cm;
  //yInit = 0.0*cm;
  //zInit =-4.8*cm;

  //with para - volume 1)
  //xInit = (-4.0*TMath::Sin(UA9Const::angleDet) + 1.1)*cm;
  //yInit = 0.0*cm;
  //zInit = -4.0*TMath::Cos(UA9Const::angleDet)*cm;

  //with para - volume 2)
  G4double lenght = 34.35/2.0*cm;
  G4double ddl = 2.0*mm;
  //xInit = (-(lenght + 0.3*cm + 0.15*cm + ddl)*TMath::Sin(UA9Const::angleDet));
  //yInit = 0.0*cm;
  //zInit = (-(lenght + 0.3*cm + 0.15*cm + ddl)*TMath::Cos(UA9Const::angleDet));

  //with para - volume into second channel
  //G4double xCopyShift = (3 + 0.5)/TMath::Cos(UA9Const::angleDet)*mm;
  //xInit = (-(5.0 - 0.3 + 0.15)*TMath::Sin(UA9Const::angleDet))*cm + 2*xCopyShift;
  //yInit = 0.0*cm;
  //zInit = (-(5.0 - 0.3 + 0.15)*TMath::Cos(UA9Const::angleDet))*cm;

  G4double LL = 4.0*cm;
  //xInit = xInit + LL;
  //yInit = yInit;
  //zInit = zInit - LL*TMath::Tan(UA9Const::angleDet);

  xInit =  10.0*cm;
  yInit =   0.0*cm;
  zInit =  36.0*cm/2.0;

  //xInit = -15.0*cm;
  //yInit = 0.0*cm;
  //zInit = 0.0*cm;

  //xInit = (-5.0*TMath::Sin(UA9Const::angleDet))*cm;
  //xInit = 0.0;
  //yInit = 1.0*cm;
  //zInit = -5.0*TMath::Cos(UA9Const::angleDet)*cm + 2.0*mm;
  //zInit = -5.0*TMath::Cos(UA9Const::angleDet)*cm + 3.0*mm;

  ///////////////////////
  _BunchXID++;
  particle = particleTable->FindParticle(_particleName);
  //particle = particleTable->FindParticle("proton");
  //G4cout<<_particleName<<G4endl;
  m = particle->GetPDGMass();
  Ekin = (TMath::Sqrt(_particleMomentum*_particleMomentum + m*m) - m);

  //with 90 deg anle - simple box
  dX = -1.00;
  dY =  0.00;
  dZ =  0.00;

  //with para - volume
  //dX =-1.0;
  //dY = 0.0;
  //dZ =-dX*TMath::Tan(UA9Const::angleDet);

  //dX = 0.0;
  //dY = -1.0;
  //dZ = 0.0;

  G4ThreeVector dir(dX, dY, dZ);
  _particleGun->SetParticleDefinition(particle);
  _particleGun->SetParticleMomentumDirection(dir);
  _particleGun->SetParticleEnergy(Ekin);  
  _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
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
