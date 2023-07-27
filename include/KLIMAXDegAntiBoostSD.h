#ifndef KLIMAXDegAntiBoostSD_h
#define KLIMAXDegAntiBoostSD_h 1

#include "SiddhartaTrackerHit.h"

#include <G4VSensitiveDetector.hh>
#include <G4SystemOfUnits.hh>

class G4Step;
class G4HCofThisEvent;

class KLIMAXDegAntiBoostSD : public G4VSensitiveDetector
{
public:
  KLIMAXDegAntiBoostSD(G4String);
  ~KLIMAXDegAntiBoostSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  SiddhartaTrackerHitsCollection* trackerCollection;
  G4double sciEnergy;
  G4double sciEnergyMax;
  G4double TimeKLDAntiBoost;
  G4double X;
  G4double Y;
  G4double Z;

  G4int kaonCounter;
  G4double kaonKinE;
};

#endif

