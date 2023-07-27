#ifndef KLIMAXCZTBoostSD_h
#define KLIMAXCZTBoostSD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>

class G4Step;
class G4HCofThisEvent;

class KLIMAXCZTBoostSD : public G4VSensitiveDetector
{
public:
  KLIMAXCZTBoostSD(G4String);
  ~KLIMAXCZTBoostSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double sciEnergy;
  G4double sciEnergyMax;
  G4double TimeKLCZTBoost;
  G4double X;
  G4double Y;
  G4double Z;

  G4int kaonCounter;
  G4int gammacounter;
  G4double kaonKinE;
};

#endif

