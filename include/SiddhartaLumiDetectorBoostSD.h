#ifndef SiddhartaLumiDetectorBoostSD_h
#define SiddhartaLumiDetectorBoostSD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>
#include <G4SystemOfUnits.hh>

class G4Step;
class G4HCofThisEvent;

class SiddhartaLumiDetectorBoostSD : public G4VSensitiveDetector
{
public:
  SiddhartaLumiDetectorBoostSD(G4String);
  ~SiddhartaLumiDetectorBoostSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double sciEnergy;
  G4double sciEnergyMax;
  G4double TimeLMBoost;
  G4double X;
  G4double Y;
  G4double Z;
};

#endif
