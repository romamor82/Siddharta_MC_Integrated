#ifndef SiddhartaAnti3SD_h
#define SiddhartaAnti3SD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>

class G4Step;
class G4HCofThisEvent;

class SiddhartaAnti3SD : public G4VSensitiveDetector
{
public:
  SiddhartaAnti3SD(G4String);
  ~SiddhartaAnti3SD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double TimeAnti[12];
  G4double EnergyAnti[12];
  G4double XYZAnti[12][3];
  G4double XYZVertexAnti[12][3];
};

#endif
