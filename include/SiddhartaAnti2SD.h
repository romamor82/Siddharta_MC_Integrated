#ifndef SiddhartaAnti2SD_h
#define SiddhartaAnti2SD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>

class G4Step;
class G4HCofThisEvent;

class SiddhartaAnti2SD : public G4VSensitiveDetector
{
public:
  SiddhartaAnti2SD(G4String);
  ~SiddhartaAnti2SD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double TimeAnti[2];
  G4double EnergyAnti[2];
  G4double XYZAnti[2][3];
  G4double XYZVertexAnti[2][3];
};

#endif
