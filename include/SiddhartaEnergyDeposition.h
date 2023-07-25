#ifndef SiddhartaEnergyDeposition_h
#define SiddhartaEnergyDeposition_h 1

#include <globals.hh>

class SiddhartaEnergyDeposition
{
public:
  SiddhartaEnergyDeposition();
  SiddhartaEnergyDeposition( const SiddhartaEnergyDeposition &right );
  SiddhartaEnergyDeposition( G4double, G4double, G4double );
  virtual ~SiddhartaEnergyDeposition();

  G4bool operator==(const SiddhartaEnergyDeposition &right) const ;
  G4bool operator< (const SiddhartaEnergyDeposition &right) const ;
  G4bool operator<=(const SiddhartaEnergyDeposition &right) const ;

  G4double GetEnergy() {return Energy;};
  G4double GetTime() {return Time;};
  G4double GetWeight() {return Weight;};

private:
  G4double Energy;
  G4double Time;
  G4double Weight;
};

#endif
