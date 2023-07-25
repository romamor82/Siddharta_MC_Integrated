#ifndef SpecialCuts_h
#define SpecialCuts_h 1

#include <G4VProcess.hh>
#include <globals.hh>
#include <G4ios.hh>

class SpecialCuts : public G4VProcess
{
public:
  SpecialCuts(const G4String& processName = "SpecialCuts");
  virtual ~SpecialCuts();

  virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track, G4double previousStepSize, G4ForceCondition* condition);
  virtual G4VParticleChange* PostStepDoIt(const G4Track& , const G4Step&);
//  no operation in  AtRestGPIL
  virtual G4double AtRestGetPhysicalInteractionLength(const G4Track& ,G4ForceCondition*) {return -1.0;};
//  no operation in  AtRestDoIt
  virtual G4VParticleChange* AtRestDoIt(const G4Track& , const G4Step& ) {return NULL;};
//  no operation in  AlongStepGPIL
  virtual G4double AlongStepGetPhysicalInteractionLength(const G4Track&, G4double , G4double , G4double& , G4GPILSelection* )
  {return -1.0;};
//  no operation in  AlongStepDoIt
  virtual G4VParticleChange* AlongStepDoIt(const G4Track& , const G4Step&) {return NULL;};

private:
  // hide assignment operator as private
     SpecialCuts& operator=(const SpecialCuts&){return *this;};
};

class MinEkineCuts : public SpecialCuts
{
public:
  MinEkineCuts(const G4String& processName = "MinEkineCuts");
  virtual ~MinEkineCuts();
// PostStep GPIL
  virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track, G4double previousStepSize, G4ForceCondition* condition);

private:
// hide assignment operator as private
  MinEkineCuts(MinEkineCuts&);
  MinEkineCuts& operator=(const MinEkineCuts& right);
};

class MaxTimeCuts : public SpecialCuts
{
public:
  MaxTimeCuts(const G4String& processName = "MaxTimeCuts");
  virtual ~MaxTimeCuts();
// PostStep GPIL
  virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track, G4double previousStepSize, G4ForceCondition* condition);

private:
// hide assignment operator as private
  MaxTimeCuts(MaxTimeCuts&);
  MaxTimeCuts& operator=(const MaxTimeCuts& right);
};

#endif
