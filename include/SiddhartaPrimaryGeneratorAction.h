#ifndef SiddhartaPrimaryGeneratorAction_h
#define SiddhartaPrimaryGeneratorAction_h 1

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4VShortLivedParticle.hh>
#include <G4SystemOfUnits.hh>

#include <iostream>
#include <fstream>

class SiddhartaDetectorConstruction;
class G4ParticleGun;
class G4Event;

const static G4double mphi_pdg = 1019.46*MeV;
const static G4double phi_Gm = 4.26*MeV;

class SiddhartaPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  SiddhartaPrimaryGeneratorAction(SiddhartaDetectorConstruction*);
  ~SiddhartaPrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*);
  
private:
  G4ParticleGun* particleGun;
  SiddhartaDetectorConstruction* myDetector;
  G4ParticleDefinition* particleVertex;

  double xx = 0;
  double xxp = 0;
  double yy = 0;
  double yyp = 0;
  double zz = 0;
  double de = 0;
  double rate = 0;
  double turn = 0;
  int done = 0;
};

#endif
