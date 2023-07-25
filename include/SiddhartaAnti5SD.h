#ifndef SiddhartaAnti5SD_h
#define SiddhartaAnti5SD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>

class G4Step;
class G4HCofThisEvent;

class SiddhartaAnti5SD : public G4VSensitiveDetector
{
public:
  SiddhartaAnti5SD(G4String);
  ~SiddhartaAnti5SD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double TimeAnti;
  G4double EnergyAnti;
  G4double XYZAnti[3];
  G4double XYZVertexAnti[3];
};

#endif
