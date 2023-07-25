#ifndef SiddhartaAnti1SD_h
#define SiddhartaAnti1SD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>

class G4Step;
class G4HCofThisEvent;

class SiddhartaAnti1SD : public G4VSensitiveDetector
{
public:
  SiddhartaAnti1SD(G4String);
  ~SiddhartaAnti1SD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double TimeAnti[2];
  G4double EnergyAnti[2];
  G4double XYZAnti[2][3];
  G4double XYZVertexAnti[2][3];

//Tables to vectors or safer structures.
// 1SD -> 5SD all the same strcutre -> inheritance instaed creation of copying the same code
};

#endif
