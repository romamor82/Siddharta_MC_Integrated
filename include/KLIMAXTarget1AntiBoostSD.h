#ifndef KLIMAXTarget1AntiBoostSD_h
#define KLIMAXTarget1AntiBoostSD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>

class G4Step;
class G4HCofThisEvent;

class KLIMAXTarget1AntiBoostSD : public G4VSensitiveDetector
{
public:
  KLIMAXTarget1AntiBoostSD(G4String);
  ~KLIMAXTarget1AntiBoostSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double sciEnergy;
  G4double sciEnergyMax;
  G4double TimeKLT1AntiBoost;
  G4double X;
  G4double Y;
  G4double Z;

  G4int kaonCounter;
  G4double kaonKinE;
};

#endif

