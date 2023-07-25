#ifndef SiddhartaKPlusSD_h
#define SiddhartaKPlusSD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>

class G4Step;
class G4HCofThisEvent;

class SiddhartaKPlusSD : public G4VSensitiveDetector
{
public:
  SiddhartaKPlusSD(G4String);
  ~SiddhartaKPlusSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double sciEnergy;
  G4double sciEnergyMax;
  G4double TimeKPlusDetector;
  G4double X;
  G4double Y;
  G4double Z;
};

#endif
