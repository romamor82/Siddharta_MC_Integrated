#ifndef SiddhartaGhostSD_h
#define SiddhartaGhostSD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>
#include <G4SystemOfUnits.hh>

class G4Step;
class G4HCofThisEvent;

class SiddhartaGhostSD : public G4VSensitiveDetector
{
public:
  SiddhartaGhostSD(G4String);
 ~SiddhartaGhostSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double TimeGhost;
  G4double X;
  G4double Y;
  G4double Z;
  G4double sciEnergy;
  G4double sciEnergyMax;
  G4double KineticGhost;
};

#endif
