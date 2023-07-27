#ifndef SiddhartaDetectorConstruction_h
#define SiddhartaDetectorConstruction_h 1

#include "SiddhartaMagneticField.h"

#include <G4VUserDetectorConstruction.hh>
#include <G4SystemOfUnits.hh>

#include <globals.hh>

class G4Tubs;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VPVParameterisation;
class G4UserLimits;
class SiddhartaDetectorMessenger;

class DegraderModificator;

class SiddhartaDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  SiddhartaDetectorConstruction();
  ~SiddhartaDetectorConstruction();

  G4VPhysicalVolume* Construct();
  G4double GetTargetFullLength() {return fTargetLength;};
  G4double GetWorldFullLength() {return fWorldLength;};

  void setTargetMaterial(G4String);
  void SetMagField(G4double);
  void SetMaxStep (G4double);

private:
  G4Box*             solidWorld;    // pointer to the solid envelope
  G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
  G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope

  G4Tubs*            solidTarget;   // pointer to the solid Target
  G4LogicalVolume*   logicTarget;   // pointer to the logical Target
  G4VPhysicalVolume* physiTarget;   // pointer to the physical Target

  G4Tubs*            solidBeamPipe;   // pointer to the solid BeamPipe
  G4LogicalVolume*   logicBeamPipe;   // pointer to the logical BeamPipe
  G4VPhysicalVolume* physiBeamPipe;   // pointer to the physical BeamPipe

  G4Box*             solidKaonDetectorTop;    // pointer to the solid KaonDetectorTop
  G4LogicalVolume*   logicKaonDetectorTop;    // pointer to the logical KaonDetectorTop
  G4VPhysicalVolume* physiKaonDetectorTop;    // pointer to the physical KaonDetectorTop

  G4Box*             solidKaonDetectorBottom;    // pointer to the solid KaonDetectorBottom
  G4LogicalVolume*   logicKaonDetectorBottom;    // pointer to the logical KaonDetectorBottom
  G4VPhysicalVolume* physiKaonDetectorBottom;    // pointer to the physical KaonDetectorBottom

  G4Box*             solidLumiDetectorBoost;    // 2020
  G4LogicalVolume*   logicLumiDetectorBoost;    // 2020
  G4VPhysicalVolume* physiLumiDetectorBoost;    // 2020

  G4Box*             solidLumiDetectorAntiboost;    // 2020
  G4LogicalVolume*   logicLumiDetectorAntiboost;    // 2020
  G4VPhysicalVolume* physiLumiDetectorAntiboost;    // 2020

  G4Box*             solidKLIMAXDegAntiBoost;    // 2020
  G4LogicalVolume*   logicKLIMAXDegAntiBoost;    // 2020
  G4VPhysicalVolume* physiKLIMAXDegAntiBoost;    // 2020

  G4Box*             solidKLIMAXTarget1AntiBoost;    // 2020
  G4LogicalVolume*   logicKLIMAXTarget1AntiBoost;    // 2020
  G4VPhysicalVolume* physiKLIMAXTarget1AntiBoost;    // 2020

  G4Box*             solidKLIMAXBoxWindowAntiBoost;    // 2020
  G4LogicalVolume*   logicKLIMAXBoxWindowAntiBoost;    // 2020
  G4VPhysicalVolume* physiKLIMAXBoxWindowAntiBoost;    // 2020

  G4Box*             solidKLIMAXTarget2AntiBoost;    // 2020
  G4LogicalVolume*   logicKLIMAXTarget2AntiBoost;    // 2020
  G4VPhysicalVolume* physiKLIMAXTarget2AntiBoost;    // 2020

  G4Box*             solidKLIMAXCZTAntiBoost;    // 2020
  G4LogicalVolume*   logicKLIMAXCZTAntiBoost;    // 2020
  G4VPhysicalVolume* physiKLIMAXCZTAntiBoost;    // 2020

  G4Box*             solidKLIMAXDegBoost;    // 2020
  G4LogicalVolume*   logicKLIMAXDegBoost;    // 2020
  G4VPhysicalVolume* physiKLIMAXDegBoost;    // 2020

  G4Box*             solidKLIMAXTargetBoost;    // 2020
  G4LogicalVolume*   logicKLIMAXTargetBoost;    // 2020
  G4VPhysicalVolume* physiKLIMAXTargetBoost;    // 2020

  G4Box*             solidKLIMAXCZTBoost;    // 2020
  G4LogicalVolume*   logicKLIMAXCZTBoost;    // 2020
  G4VPhysicalVolume* physiKLIMAXCZTBoost;    // 2020

  G4Box*             solidKPlusDetector;    // pointer to the solid KaonPlusDetector
  G4LogicalVolume*   logicKPlusDetector;    // pointer to the logical KaonPlusDetector
  G4VPhysicalVolume* physiKPlusDetector;    // pointer to the physical KaonPlusDetector

  G4Box*             solidDegrader;    // pointer to the solid KaonPlusDetector
  G4LogicalVolume*   logicDegrader;    // pointer to the logical KaonPlusDetector
  G4VPhysicalVolume* physiDegrader;    // pointer to the physical KaonPlusDetector

  G4Tubs*            solidShielding;    // pointer to the solid KaonPlusDetector
  G4LogicalVolume*   logicShielding;    // pointer to the logical KaonPlusDetector
  G4VPhysicalVolume* physiShielding;    // pointer to the physical KaonPlusDetector

  G4Material*         TargetMater;  // pointer to the target material
  G4Material*         BeamPipeMater;  // pointer to the BeamPipe material
  G4Material*         KaonDetectorTopMater;  // pointer to the KaonDetectorTop material
  G4Material*         KaonDetectorBottomMater;  // pointer to the KaonDetectorBottom material
  G4Material*         KPlusDetectorMater;  // pointer to the KPlusDetector material
  G4Material*         DegraderMater;  // pointer to the KPlusDetector material
  G4Material*         ShieldingMater;  // pointer to the KPlusDetector material

  G4UserLimits* stepLimit;             // pointer to user step limits

  SiddhartaMagneticField* fpMagField;   // pointer to the magnetic field

  SiddhartaDetectorMessenger* detectorMessenger;  // pointer to the Messenger

  G4double fWorldLength;            // Full length of the world volume
  G4double fTargetLength;           // Full length of Target

  DegraderModificator* fDegModif;
};

class DegraderModificator
{
public:
  DegraderModificator() {
    fStepMaxThickness = {0., 0., 1000.*um, 100.*um, 100.*um, 150.*um, 100.*um, 100.*um, 100.*um, 100.*um, 0., 0.};
    fDegThicknessFromHalf = 550.*um;
  };
  void SetThicknessFromHaLF(G4double thick) {fDegThicknessFromHalf = thick;};
  void SetMaxThicknessAtLayer(G4double thick, int id) {
    if (id >= 1 && id <= fStepMaxThickness.size()) {
      fStepMaxThickness.at(id-1) = thick;
    } else {
      G4cout << "DegraderModificator error: SetThicknessAtLayer. Trying to set " << id;
      G4cout << " layer, but there are only layers from 1 to " << fStepMaxThickness.size() << G4endl;
    }
  };
  
  G4double GetThickAtLayer(int id) {
    if (id >= 1 && id <= fStepMaxThickness.size()) {
      std::vector<G4double> fStepThickness = fStepMaxThickness;
      
      if (fDegThicknessFromHalf > 0) {
	if (id > 0.5*fStepMaxThickness.size()) {
	  return (fStepMaxThickness.at(id-1) == 0 ? 0.1*um : fStepMaxThickness.at(id-1));
	}

	G4double remainingThickAtHalf = fDegThicknessFromHalf;
	int currentID = (int)(0.5*fStepMaxThickness.size());
	G4double resultAtID = 0;
	bool whileTest = true;
	while (whileTest && currentID >= 0) {
	  if (remainingThickAtHalf < fStepMaxThickness.at(currentID-1)) {
 	    if (id != currentID) {
	      resultAtID = 0;
	    } else {
	      resultAtID = remainingThickAtHalf;
  	    }
	    whileTest = false;
	  } else {
 	    if (id == currentID) {
	      resultAtID = fStepMaxThickness.at(currentID-1);
	      whileTest = false;
	    }
          }
 	  if (resultAtID < 1e-10)
            resultAtID = 0.1*um;
	  remainingThickAtHalf -= fStepMaxThickness.at(currentID-1);
	  currentID--;
	}
	return resultAtID;
      } else {
	if (id <= 0.5*fStepMaxThickness.size()) {
	  return 0.1*um;
	}

	G4double remainingThickAtHalf = -fDegThicknessFromHalf;
	int currentID = (int)(0.5*fStepMaxThickness.size()) + 1;
	G4double resultAtID = 0;
	bool whileTest = true;

	while (whileTest && currentID <= fStepMaxThickness.size()) {
	  if (remainingThickAtHalf < fStepMaxThickness.at(currentID-1)) {
 	    if (id != currentID) {
	      resultAtID = fStepMaxThickness.at(id-1);
	    } else {
	      resultAtID = fStepMaxThickness.at(id-1) - remainingThickAtHalf;
  	    }
	    whileTest = false;
	  } else {
 	    if (id == currentID) {
	      resultAtID = 0;
	      whileTest = false;
	    }
          }
 	  if (resultAtID < 1e-10)
            resultAtID = 0.1*um;
	  remainingThickAtHalf -= fStepMaxThickness.at(currentID-1);
	  currentID++;
	}
	return resultAtID;
      }
    } else {
      G4cout << "DegraderModificator error: GetThickAtLayer. Trying to get " << id;
      G4cout << " layer, but there are only layers from 1 to " << fStepMaxThickness.size() << G4endl;
    }
    return 0.1*um;
  };

  int GetThickFromHalf() {return fDegThicknessFromHalf;};
  int GetLayerSize() {return fStepMaxThickness.size();};

private:
  std::vector<G4double> fStepMaxThickness;
  G4double fDegThicknessFromHalf;
};

#endif
