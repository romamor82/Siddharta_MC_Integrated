#include "../include/SpecialCuts.h"

#include <G4VParticleChange.hh>
#include <G4Track.hh>
#include <G4Step.hh>

SpecialCuts::SpecialCuts(const G4String& aName) : G4VProcess(aName)
{
  if (verboseLevel>1) {
    G4cout << GetProcessName() << " is created "<< G4endl;
  }
}

SpecialCuts::~SpecialCuts()
{
}

G4VParticleChange* SpecialCuts::PostStepDoIt(const G4Track& aTrack, const G4Step& )
//
// Stop the current particle, if requested by G4UserLimits
//
{
  aParticleChange.Initialize(aTrack);
  aParticleChange.ProposeEnergy(0.) ;
  aParticleChange.ProposeLocalEnergyDeposit (aTrack.GetKineticEnergy()) ;
  aParticleChange.ProposeTrackStatus(fStopButAlive);
  return &aParticleChange;
}

G4double SpecialCuts::PostStepGetPhysicalInteractionLength(const G4Track& , G4double , G4ForceCondition* )
{
  return DBL_MAX;
}

#include <G4PhysicalConstants.hh>
#include <G4LossTableManager.hh>
#include <G4UserLimits.hh>

MinEkineCuts::MinEkineCuts(const G4String& aName) : SpecialCuts(aName)
{
  if (verboseLevel>1) {
    G4cout << GetProcessName() << " is created "<< G4endl;
  }
  SetProcessType(fUserDefined);
}

MinEkineCuts::~MinEkineCuts()
{}

MinEkineCuts::MinEkineCuts(MinEkineCuts&) : SpecialCuts()
{}

G4double MinEkineCuts::PostStepGetPhysicalInteractionLength(const G4Track& aTrack, G4double , G4ForceCondition* condition)
{
// condition is set to "Not Forced"
  *condition = NotForced;

  G4double proposedStep = DBL_MAX;
// get the pointer to UserLimits
  G4UserLimits* pUserLimits = aTrack.GetVolume()->GetLogicalVolume()->GetUserLimits();
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4ParticleDefinition* aParticleDef = aTrack.GetDefinition();

  if (pUserLimits && aParticleDef->GetPDGCharge() != 0.0) {
//min kinetic energy
    G4double temp = DBL_MAX;
    G4double eKine = aParticle->GetKineticEnergy();
    const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
    G4double eMin = pUserLimits->GetUserMinEkine(aTrack);

    if (eKine < eMin)
      return 0.;

    G4double rangeNow = DBL_MAX;
    G4LossTableManager* lossManager = G4LossTableManager::Instance();
    rangeNow = lossManager->GetRange(aParticleDef, eKine, couple);
// charged particles only
    G4double rangeMin = lossManager->GetRange(aParticleDef, eMin, couple);
    temp = rangeNow - rangeMin;

    if (proposedStep > temp)
      proposedStep = temp;
  }
  return proposedStep;
}

MaxTimeCuts::MaxTimeCuts(const G4String& aName) : SpecialCuts(aName)
{
  if (verboseLevel>1) {
    G4cout << GetProcessName() << " is created "<< G4endl;
  }
  SetProcessType(fUserDefined);
}

MaxTimeCuts::~MaxTimeCuts()
{}

MaxTimeCuts::MaxTimeCuts(MaxTimeCuts&) : SpecialCuts()
{}


G4double MaxTimeCuts::PostStepGetPhysicalInteractionLength(const G4Track& aTrack, G4double , G4ForceCondition* condition)
{
// condition is set to "Not Forced"
  *condition = NotForced;

   G4double proposedStep = DBL_MAX;
// get the pointer to UserLimits
   G4UserLimits* pUserLimits = aTrack.GetVolume()->GetLogicalVolume()->GetUserLimits();
   const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

   if (pUserLimits) {
     G4double temp = DBL_MAX;
//max time limit
     G4double dTime = (pUserLimits->GetUserMaxTime(aTrack) - aTrack.GetGlobalTime());
     if (dTime < 0. ) {
       proposedStep = 0.;
     } else {
       G4double beta = (aParticle->GetTotalMomentum())/(aParticle->GetTotalEnergy());
       temp = beta*c_light*dTime;

       if (proposedStep > temp)
        proposedStep = temp;
     }
   }
   return proposedStep;
}
