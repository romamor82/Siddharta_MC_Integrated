#ifndef KLIMAXBoxWindowAntiBoostSD_h
#define KLIMAXBoxWindowAntiBoostSD_h 1

#include "G4VSensitiveDetector.hh"
#include "SiddhartaTrackerHit.h"

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class KLIMAXBoxWindowAntiBoostSD : public G4VSensitiveDetector
{
  public:
      KLIMAXBoxWindowAntiBoostSD(G4String);
     ~KLIMAXBoxWindowAntiBoostSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);

  private:
      SiddhartaTrackerHitsCollection* trackerCollection;
      G4double sciEnergy;
      G4double sciEnergyMax;
      G4double TimeKLBWAntiBoost;
      G4double X;
      G4double Y;
      G4double Z;

      G4int kaonCounter; 
      G4double kaonKinE; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

