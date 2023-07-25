#include "../include/SiddhartaMagneticField.h"

#include <G4TransportationManager.hh>
#include <G4FieldManager.hh>

SiddhartaMagneticField::SiddhartaMagneticField() : G4UniformMagField(G4ThreeVector())
{
  GetGlobalFieldManager()->SetDetectorField(this);
  GetGlobalFieldManager()->CreateChordFinder(this);
}

SiddhartaMagneticField::SiddhartaMagneticField(G4ThreeVector fieldVector) : G4UniformMagField(fieldVector)
{
  GetGlobalFieldManager()->SetDetectorField(this);
  GetGlobalFieldManager()->CreateChordFinder(this);
}

SiddhartaMagneticField::~SiddhartaMagneticField() {}

void SiddhartaMagneticField::SetMagFieldValue(G4double fieldValue)
{
  SetMagFieldValue(G4ThreeVector(fieldValue, 0, 0));
}

void SiddhartaMagneticField::SetMagFieldValue(G4ThreeVector fieldVector)
{
  G4FieldManager* fieldMgr = GetGlobalFieldManager();

  if (fieldVector != G4ThreeVector(0.,0.,0.)) {
    SetFieldValue(fieldVector);
    fieldMgr->SetDetectorField(this);
  } else {
    G4MagneticField* magField = 0;
    fieldMgr->SetDetectorField(magField);
  }
}

G4FieldManager* SiddhartaMagneticField::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()->GetFieldManager();
}
