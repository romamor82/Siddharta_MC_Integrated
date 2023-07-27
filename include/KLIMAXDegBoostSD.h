#ifndef KLIMAXDegBoostSD_h
#define KLIMAXDegBoostSD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>
#include <G4SystemOfUnits.hh>

class G4Step;
class G4HCofThisEvent;

class KLIMAXDegBoostSD : public G4VSensitiveDetector
{
public:
  KLIMAXDegBoostSD(G4String);
  ~KLIMAXDegBoostSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double sciEnergy;
  G4double sciEnergyMax;
  G4double TimeKLDBoost;
  G4double X;
  G4double Y;
  G4double Z;

  G4int kaonCounter;
  G4double kaonKinE;
};

#endif

