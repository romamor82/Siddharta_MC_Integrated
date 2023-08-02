#ifndef KLIMAXHPGeBoostSD_h
#define KLIMAXHPGeBoostSD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>

class G4Step;
class G4HCofThisEvent;

class KLIMAXHPGeBoostSD : public G4VSensitiveDetector
{
public:
  KLIMAXHPGeBoostSD(G4String);
  ~KLIMAXHPGeBoostSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double sciEnergy;
  G4double sciEnergyMax;
  G4double TimeKLHPGeBoost;
  G4double X;
  G4double Y;
  G4double Z;

  G4int kaonCounter;
  G4int gammacounter;
  G4double kaonKinE;
};

#endif

