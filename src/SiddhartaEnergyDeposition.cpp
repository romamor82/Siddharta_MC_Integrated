#include "../include/SiddhartaEnergyDeposition.h"

SiddhartaEnergyDeposition::SiddhartaEnergyDeposition() {}

SiddhartaEnergyDeposition::SiddhartaEnergyDeposition(G4double energy, G4double time, G4double weight) : Energy(energy), Time(time), Weight(weight)
{}

SiddhartaEnergyDeposition::SiddhartaEnergyDeposition(const SiddhartaEnergyDeposition &right) : Energy(right.Energy), Time(right.Time), Weight(right.Weight)
{}

SiddhartaEnergyDeposition::~SiddhartaEnergyDeposition() {}

G4bool SiddhartaEnergyDeposition::operator==(const SiddhartaEnergyDeposition &right) const
{
  return Time == right.Time;
}

G4bool SiddhartaEnergyDeposition::operator<(const SiddhartaEnergyDeposition &right) const
{
  return Time < right.Time;
}

G4bool SiddhartaEnergyDeposition::operator<=(const SiddhartaEnergyDeposition &right) const
{
  return Time <= right.Time;
}

