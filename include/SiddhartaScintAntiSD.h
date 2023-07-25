#ifndef SiddhartaScintAntiSD_h
#define SiddhartaScintAntiSD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>
#include <G4SystemOfUnits.hh>

class G4Step;
class G4HCofThisEvent;

class SiddhartaScintAntiSD : public G4VSensitiveDetector
{
public:
  SiddhartaScintAntiSD(G4String);
  ~SiddhartaScintAntiSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4int nHit[48];
  G4int nHitFull;
  G4double TimeScintAnti[700];
  G4double EnergyDep[700];
  G4int Copy[700];
};

#endif
