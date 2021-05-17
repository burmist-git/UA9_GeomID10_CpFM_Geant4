#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4String.hh"

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);
  
public:
  void SetParticleName(G4String particleName) {_particleName = particleName;}
  G4String GetParticleName() {return _particleName;}
  void SetParticleMomentum(G4double momentum) {_particleMomentum = momentum;}
  G4double GetParticleMomentum() {return _particleMomentum;}
  void SetPhiAngle(G4double val) {_PhiAngle = val;}
  G4double GetPhiAngle() { return _PhiAngle; }
  void SetThetaAngle(G4double val) {_ThetaAngle = val;}
  G4double GetThetaAngle() { return _ThetaAngle; }
  void SetSinglePhoton(G4bool singlePhoton) { _singlePhoton = singlePhoton;}
  G4bool SinglePhotonGenerator() { return _singlePhoton;}
  G4int GetBunchXID(){return _BunchXID;};
private:
  G4ParticleGun* _particleGun;
  G4String _particleName;
  G4double _particleMomentum;
  G4double _PhiAngle;
  G4double _ThetaAngle;

  G4bool _singlePhoton;

  G4int _BunchXID;

  G4int GenFlatInt(G4int iMin,G4int iMax);
  void generateThetaAndPhi();
  G4double genCos2dist();

};

#endif
