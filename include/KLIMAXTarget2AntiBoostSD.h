#ifndef KLIMAXTarget2AntiBoostSD_h
#define KLIMAXTarget2AntiBoostSD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>

class G4Step;
class G4HCofThisEvent;

class KLIMAXTarget2AntiBoostSD : public G4VSensitiveDetector
{
public:
  KLIMAXTarget2AntiBoostSD(G4String);
  ~KLIMAXTarget2AntiBoostSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double sciEnergy;
  G4double sciEnergyMax;
  G4double TimeKLT2AntiBoost;
  G4double X;
  G4double Y;
  G4double Z;

  G4int kaonCounter;
  G4double kaonKinE;
};

#endif

