#ifndef SiddhartaLumiDetectorAntiboostSD_h
#define SiddhartaLumiDetectorAntiboostSD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>
#include <G4SystemOfUnits.hh>

class G4Step;
class G4HCofThisEvent;

class SiddhartaLumiDetectorAntiboostSD : public G4VSensitiveDetector
{
public:
  SiddhartaLumiDetectorAntiboostSD(G4String);
  ~SiddhartaLumiDetectorAntiboostSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double sciEnergy;
  G4double sciEnergyMax;
  G4double TimeLMAntiboost;
  G4double X;
  G4double Y;
  G4double Z;

  G4int kaonCounter;
};

#endif
