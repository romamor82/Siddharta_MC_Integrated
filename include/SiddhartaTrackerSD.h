#ifndef SiddhartaTrackerSD_h
#define SiddhartaTrackerSD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>
#include <G4SystemOfUnits.hh>

#include <map>

class G4Step;
class G4HCofThisEvent;

class SiddhartaTrackerSD : public G4VSensitiveDetector
{
public:
  SiddhartaTrackerSD(G4String);
  ~SiddhartaTrackerSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4int nHit[48][8];
  G4int track[48][8];
  G4double sddEnergy[48][8];
  G4int nHitFull;
  G4int sddNb[48*8];
  G4int chipNb[48*8];
  std::map<G4int,G4int> trackhit;
};

#endif
