//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

#include "../include/G4KaonMinusAbsorptionAtRest.h"
#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaHisto.h"

//#include <G4StopDeexcitationAlgorithm.hh>
//#include <G4StopTheoDeexcitation.hh>
//#include <G4StopDeexcitation.hh>
//
//equivalent to G4PhaseSpaceDecayChannel
//#include <G4ReactionKinematics.hh>

#include <G4HadronicProcessStore.hh>

#define G4RandBreitWigner CLHEP::RandBreitWigner

G4KaonMinusAbsorptionAtRest::G4KaonMinusAbsorptionAtRest(const G4String& processName, G4ProcessType aType) : G4VRestProcess (processName, aType)
{
  if (verboseLevel > 0) {
    G4cout << GetProcessName() << " is created "<< G4endl;
  }
  SetProcessSubType(fHadronAtRest);

  // see Cohn et al, PLB27(1968) 527;
  //     Davis et al, PLB1(1967) 434;
  pionAbsorptionRate = 0.07;

  // see  VanderVelde-Wilquet et al, Nuov.Cim.39A(1978)538;
  // see  VanderVelde-Wilquet et al, Nuov.Cim.38A(1977)178;
  // see  VanderVelde-Wilquet et al, Nucl.Phys.A241(1975)511;
  // primary production rates ( for absorption on Carbon)
  // .. other elements are extrapolated by the halo factor.
  rateLambdaZeroPiZero = 0.052;
  rateSigmaMinusPiPlus = 0.199;
  rateSigmaPlusPiMinus = 0.446;
  rateSigmaZeroPiZero  = 0.303;
  rateLambdaZeroPiMinus = 0.568;
  rateSigmaZeroPiMinus  = 0.216;
  rateSigmaMinusPiZero  = 0.216;

  // for sigma- p -> lambda n
  //     sigma+ n -> lambda p
  //     sigma- n -> lambda
  // all values compatible with 0.55 same literature as above.
  sigmaPlusLambdaConversionRate = 0.55;
  sigmaMinusLambdaConversionRate = 0.55;
  sigmaZeroLambdaConversionRate = 0.55;

  //atomic cascade transitions
  //see yields from DEAR-SIDDHARTA bibliography
  //M.Iliescu, Feb. 2010
  yieldKa=1.;probKa=.3;
  yieldKb=1.;probKb=.2;
  yieldKg=1.;probKg=.15;
  yieldKd=1.;probKd=.1;
  yieldKe=1.;probKe=.1;
  yieldKz=1.;probKz=.1;
  yieldKi=1.;probKi=.05;

  yieldXaOther=1.;probXaOther=.75;
  G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);
}


G4KaonMinusAbsorptionAtRest::~G4KaonMinusAbsorptionAtRest()
{
  G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
}

void G4KaonMinusAbsorptionAtRest::PreparePhysicsTable(const G4ParticleDefinition& p)
{
  G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this, &p);
}

void G4KaonMinusAbsorptionAtRest::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

G4VParticleChange* G4KaonMinusAbsorptionAtRest::AtRestDoIt(const G4Track& track, const G4Step& )
{
  stoppedHadron = track.GetDynamicParticle();

  if (!IsApplicable(*(stoppedHadron->GetDefinition()))) {
    G4cerr << "G4KaonMinusAbsorptionAtRest:ERROR, particle must be a Kaon!" << G4endl;
    return nullptr;
  }
  currentMaterial = track.GetMaterial();
  G4double currentTime = track.GetGlobalTime();

  nucleus = nullptr;
  do {
    nucleus = new G4Nucleus(currentMaterial);

// cascade modified, original value 1.5
    if (nucleus->GetA_asInt() < 0.5) {
      delete nucleus;
      nucleus = nullptr;
    }

  } while (nucleus == nullptr);

  G4int Z = nucleus->GetZ_asInt();
  G4int A = nucleus->GetA_asInt();
  G4DynamicParticleVector* absorptionProducts = KaonNucleonReaction();

// For A=3 N=0 fragmentation case in He target 12.2020 H. Shi
  G4int fragZ = nucleus->GetZ_asInt();
  G4int fragA = nucleus->GetA_asInt();
  if (A >= 1.5 && !(fragA == 3 && fragZ == 0) )     // not removing the A=3 Z=0 error case..
  {
//Secondary interactions
    G4DynamicParticle* thePion;
    for (int i=0; i<absorptionProducts->size(); i++) {
      thePion = (*absorptionProducts)[i];
      if (thePion->GetDefinition() == G4PionMinus::PionMinus()
            || thePion->GetDefinition() == G4PionPlus::PionPlus()
            || thePion->GetDefinition() == G4PionZero::PionZero())
      {
        if (AbsorbPionByNucleus(thePion)) {
          absorptionProducts->erase(absorptionProducts->begin()+i);
          i--;
          delete thePion;
          if (verboseLevel > 1)
            G4cout << "G4KaonMinusAbsorption::AtRestDoIt: Pion absorbed in Nucleus" << G4endl;
        }
      }
    }

    G4DynamicParticle* theSigma;
    G4DynamicParticle* theLambda;
    for (int i=0; i<absorptionProducts->size(); i++) {
      theSigma = (*absorptionProducts)[i];
      if (theSigma->GetDefinition() == G4SigmaMinus::SigmaMinus()
            || theSigma->GetDefinition() == G4SigmaPlus::SigmaPlus()
            || theSigma->GetDefinition() == G4SigmaZero::SigmaZero())
      {
        theLambda = SigmaLambdaConversion(theSigma);
        if (theLambda != 0) {
          absorptionProducts->erase(absorptionProducts->begin()+i);
          i--;
          delete theSigma;
          absorptionProducts->push_back(theLambda);

          if (verboseLevel > 1)
            G4cout << "G4KaonMinusAbsorption::AtRestDoIt: SigmaLambdaConversion Done" << G4endl;
        }
      }
    }

// Nucleus deexcitation
    G4double productEnergy = 0.;
    G4ThreeVector pProducts(0.,0.,0.);

    unsigned int nAbsorptionProducts = 0;
    if (absorptionProducts != nullptr)
      nAbsorptionProducts = absorptionProducts->size();

    for (int i=0; i<nAbsorptionProducts; i++) {
      pProducts += (*absorptionProducts)[i]->GetMomentum();
      productEnergy += (*absorptionProducts)[i]->GetKineticEnergy();
    }

    G4int newZ = nucleus->GetZ_asInt();
    G4int newA = nucleus->GetA_asInt();

    G4double bDiff = G4NucleiProperties::GetBindingEnergy(A,Z) -
                            G4NucleiProperties::GetBindingEnergy(newA, newZ);

    G4StopDeexcitationAlgorithm* nucleusAlgorithm = new G4StopTheoDeexcitation();
    G4StopDeexcitation stopDeexcitation(nucleusAlgorithm);

    nucleus->AddExcitationEnergy(bDiff);

    G4double energyDeposit = nucleus->GetEnergyDeposit();
    if (verboseLevel > 0) {
      G4cout << " -- KaonAtRest -- excitation = " << energyDeposit
             << ", pNucleus = " << pProducts
             << ", A: " << A
             << ", " << newA
             << ", Z: " << Z
             << ", " << newZ << G4endl;
    }

    if (energyDeposit < 0.) {
      G4Exception("G4KaonMinusAbsorptionAtRest", "007", FatalException, "AtRestDoIt -- excitation energy < 0");
    }
    delete nucleus;

    if (newA == 3 && newZ == 0) {
      G4cout << "---- Abnormal fragmentation products A = 3  Z = 0 -----" << G4endl;
      G4cout << "     Old A = " << A << "  old Z = " << Z << G4endl;
      G4cout << "     NewA = "  << newA << "  newZ = " << newZ << G4endl;
      G4cout << "     Skip stopDeexcitation break up     " << G4endl;
      G4cout << "     After KaonNucleonReaction(): A = "  << fragA << "  Z = " << fragZ << G4endl;
    } else {
      G4ReactionProductVector* fragmentationProducts = stopDeexcitation.DoBreakUp(newA, newZ, energyDeposit, pProducts);

      unsigned nFragmentationProducts = 0;
      if (fragmentationProducts != 0)
        nFragmentationProducts = fragmentationProducts->size();

//Initialize ParticleChange -> Internal variable in G4ParticleChange <- Note in Geant4 notes as possible to remove!!!
      aParticleChange.Initialize(track);
      aParticleChange.SetNumberOfSecondaries(G4int(nAbsorptionProducts + nFragmentationProducts));

//Update List of alive particles. put energy deposit at the right place ...
      for (int i=0; i<nAbsorptionProducts; i++) {
        aParticleChange.AddSecondary((*absorptionProducts)[i], currentTime);
      }
      if (absorptionProducts != nullptr)
        delete absorptionProducts;

      for(int i=0; i<nFragmentationProducts; i++) {
        G4DynamicParticle * aNew = new G4DynamicParticle((*fragmentationProducts)[i]->GetDefinition(),
                                                            (*fragmentationProducts)[i]->GetTotalEnergy(),
                                                            (*fragmentationProducts)[i]->GetMomentum());
        G4double newTime = aParticleChange.GetGlobalTime((*fragmentationProducts)[i]->GetFormationTime());
        aParticleChange.AddSecondary(aNew, newTime);
        delete (*fragmentationProducts)[i];
      }
      if (fragmentationProducts != 0)
        delete fragmentationProducts;
    }
  } else { //else works for -> A<1.5 || fragA==3 && fragZ==0
    unsigned nAbsorptionProducts = 0;

    delete nucleus;
    if (absorptionProducts != nullptr)
      nAbsorptionProducts = absorptionProducts->size();

//Initialize ParticleChange -> Internal variable in G4ParticleChange <- Note in Geant4 notes as possible to remove!!!
    aParticleChange.Initialize(track);
    aParticleChange.SetNumberOfSecondaries(G4int(nAbsorptionProducts) );

//Update List of alive particles.
    for (int i=0; i<nAbsorptionProducts; i++) {
      aParticleChange.AddSecondary((*absorptionProducts)[i],currentTime);
    }
    if (absorptionProducts != nullptr)
      delete absorptionProducts;
  }
  aParticleChange.ProposeTrackStatus(fStopAndKill); //fStopAndKill -> internal state (enum) in G4TrackStatus
  return &aParticleChange;
}


G4DynamicParticle G4KaonMinusAbsorptionAtRest::GetAbsorbingNucleon()
{
  G4DynamicParticle aNucleon;
//Atomic cascade modified
  if (nucleus->GetA_asInt() >= 1.5) {
//Get nucleon definition, based on Z,N of current Nucleus
  	aNucleon.SetDefinition(SelectAbsorbingNucleon());

//Fermi momentum distribution in three dimensions
  	G4ThreeVector pFermi = nucleus->GetFermiMomentum();
  	aNucleon.SetMomentum(pFermi);
  } else {
    aNucleon.SetDefinition(G4Proton::Proton());  //kaonic hydrogen case
    aNucleon.SetMomentum(G4ThreeVector(0.,0.,0.));
  }
  return aNucleon;
}

G4ParticleDefinition* G4KaonMinusAbsorptionAtRest::SelectAbsorbingNucleon()
{
// (Ch. Voelcker) extended from ReturnTargetParticle():
// Choose a proton or a neutron as the absorbing particle,
// taking weight into account!
// Update nucleon's atomic numbers.
  G4ParticleDefinition* absorbingParticleDef;
  G4double ranflat = G4UniformRand();

  G4int myZ = nucleus->GetZ_asInt();   // number of protons
  G4int myN = nucleus->GetA_asInt();   // number of nucleons (not neutrons!!)

// See  VanderVelde-Wilquet et al, Nuov.Cim.39A(1978)538;
  G4double carbonRatioNP = 0.18;  // (Rn/Rp)c, see page 544

  G4double neutronProtonRatio = NeutronHaloFactor(myZ,myN)*carbonRatioNP*(double)(myN-myZ)/(double)myZ;
  G4double protonProbability = 1./(1. + neutronProtonRatio);

  if (ranflat < protonProbability) {
    absorbingParticleDef = G4Proton::Proton();
    myZ -= 1.;
  } else {
    absorbingParticleDef = G4Neutron::Neutron();
  }
//It may call an exception. Is it inteded?
  myN -= 1.;
// remove the interacting nucleon from the current nucleus
  if (myN < 1) {
    myN = 1;
  }
  if (myZ < 0) {
    myZ = 0;
  }
  nucleus->SetParameters(myN,myZ);
  return absorbingParticleDef;
}


G4double G4KaonMinusAbsorptionAtRest::NeutronHaloFactor(G4double Z, G4double N)
{
// this function should take care of the probability for absorption
// on neutrons, depending on number of protons Z and number of neutrons N-Z
// parametrisation from fit to
// VanderVelde-Wilquet et al, Nuov.Cim.39A(1978)538;

  if (Z == 1.)
    return 1.389; // deuterium
  else if (Z == 2.)
    return 1.78; // helium
  else if (Z == 10.)
    return 0.66; // neon
  else
    return 0.6742+(N-Z)*0.06524;
}


G4DynamicParticleVector* G4KaonMinusAbsorptionAtRest::KaonNucleonReaction()
{
  G4DynamicParticleVector* products = new G4DynamicParticleVector();

  G4double ranflat = G4UniformRand();
  G4double prob = 0;

  G4ParticleDefinition* producedBaryonDef;
  G4ParticleDefinition* producedMesonDef;
  G4ParticleDefinition* producedBosonDef;  // atomic cascade photon

  G4int iniZ = nucleus->GetZ_asInt();
  G4int iniA = nucleus->GetA_asInt();

  G4DynamicParticle aNucleon = GetAbsorbingNucleon();
  G4double nucleonMass;

  if (aNucleon.GetDefinition() == G4Proton::Proton()) {
    nucleonMass = proton_mass_c2 + electron_mass_c2;
    if ((prob += rateLambdaZeroPiZero) > ranflat) {
      producedBaryonDef = G4Lambda::Lambda();
      producedMesonDef  = G4PionZero::PionZero();
    } else if ((prob += rateSigmaPlusPiMinus) > ranflat) {
      producedBaryonDef = G4SigmaPlus::SigmaPlus();
      producedMesonDef  = G4PionMinus::PionMinus();
    } else if ((prob += rateSigmaMinusPiPlus) > ranflat) {
      producedBaryonDef = G4SigmaMinus::SigmaMinus();
      producedMesonDef  = G4PionPlus::PionPlus();
    } else {
      producedBaryonDef = G4SigmaZero::SigmaZero();
      producedMesonDef  = G4PionZero::PionZero();
    }
  } else if (aNucleon.GetDefinition() == G4Neutron::Neutron()) {
    nucleonMass = neutron_mass_c2;
    if ((prob += rateLambdaZeroPiMinus) > ranflat) {
      producedBaryonDef = G4Lambda::Lambda();
      producedMesonDef  = G4PionMinus::PionMinus();
    } else if ((prob += rateSigmaZeroPiMinus) > ranflat) {
      producedBaryonDef = G4SigmaZero::SigmaZero();
      producedMesonDef = G4PionMinus::PionMinus();
    } else {
      producedBaryonDef = G4SigmaMinus::SigmaMinus();
      producedMesonDef  = G4PionZero::PionZero();
    }
  } else {
    if (verboseLevel > 0) {
      G4cout << "G4KaonMinusAbsorption::KaonNucleonReaction: " << aNucleon.GetDefinition()->GetParticleName()
             << " is not a good nucleon - check G4Nucleus::ReturnTargetParticle()!" << G4endl;
    }
    return nullptr;
  }
  G4DynamicParticle modifiedHadron = (*stoppedHadron);

  if (iniA >= 1.5) {            // kaonic hydrogen cascade
    G4int newZ = nucleus->GetZ_asInt();
    G4int newA = nucleus->GetA_asInt();

// Modify the Kaon mass to take nuclear binding energy into account
// .. using mass formula ..
// .. using mass table ..
// equivalent to '-initialBindingEnergy+nucleus.GetBindingEnergy' !

    G4double nucleonBindingEnergy = G4NucleiProperties::GetBindingEnergy(newA, newZ) - G4NucleiProperties::GetBindingEnergy(iniA, iniZ);
    modifiedHadron.SetMass(stoppedHadron->GetMass() + nucleonBindingEnergy);
  }
//Is it needed still?
//else {							//implicit
//	modifiedHadron.SetMass(stoppedHadron->GetMass());
//}

// Setup outgoing dynamic particles
  G4ThreeVector dummy(0.,0.,0.);
  G4DynamicParticle* producedBaryon = new G4DynamicParticle(producedBaryonDef, dummy);
  G4DynamicParticle* producedMeson = new G4DynamicParticle(producedMesonDef, dummy);

// Produce the secondary particles in a twobody process:
  G4ReactionKinematics theReactionKinematics;
  theReactionKinematics.TwoBodyScattering(&modifiedHadron, &aNucleon, producedBaryon, producedMeson);

//====kaonic atom cascade photons=================================================


  // ******************  OLD *****************************//

/*  
  producedBosonDef = G4Gamma::Gamma();
  G4double photonEnergy;
  unsigned nphotons = 0;

  std::vector<G4double> photonEnergy1;
  if (((iniA) == 27) && ((iniZ) == 13)) {               //Aluminium
    photonEnergy1.push_back(10435.1*eV);
    photonEnergy1.push_back(17840.0*eV);
    photonEnergy1.push_back(5113.4*eV);
    photonEnergy1.push_back(9025.0*eV);

// added 16.06.2021 H. Shi KAl 7-6 from KG eq. calc.
    photonEnergy1.push_back(16058.0*eV);
    photonEnergy1.push_back(26518.4*eV); // added Jan 2023
    photonEnergy1.push_back(7150.9*eV);
    photonEnergy1.push_back(12262.8*eV);
  } else if (((iniA) == 48) && ((iniZ) == 22)) {        //Titanium
    photonEnergy1.push_back(8312.0*eV);
    photonEnergy1.push_back(14770.0*eV);
    photonEnergy1.push_back(6467.0*eV);
    photonEnergy1.push_back(11590.0*eV);
  } else if (((iniA) == 12) && ((iniZ) == 6)) {         //Carbon
    photonEnergy1.push_back(10216.5*eV);
    photonEnergy1.push_back(15809.0*eV);
    photonEnergy1.push_back(5544.9*eV);
    photonEnergy1.push_back(8885.8*eV);
  } else if (((iniA) == 16) && ((iniZ) == 8)) {         //Oxygen
    photonEnergy1.push_back(9968.7*eV);
    photonEnergy1.push_back(16062.0*eV);
    photonEnergy1.push_back(6006.8*eV);
    photonEnergy1.push_back(9958.0*eV);
  } else if (((iniA) == 14) && ((iniZ) == 7)) {         //Nitrogen
    photonEnergy1.push_back(13995.4*eV);
    photonEnergy1.push_back(21588.2*eV);
    photonEnergy1.push_back(7595.4*eV);
    photonEnergy1.push_back(12223*eV);
    photonEnergy1.push_back(4577.1*eV);
    photonEnergy1.push_back(7578*eV);
  } else if (((iniA) == 19) && ((iniZ) == 9)) {         //Fluorine
    photonEnergy1.push_back(7642.7*eV);
    photonEnergy1.push_back(12599*eV);
  } else if (((iniA) == 28) && ((iniZ) == 14)) {        //Silicon
    photonEnergy1.push_back(8300.1*eV);
    photonEnergy1.push_back(14233.1*eV);
    photonEnergy1.push_back(5935.0*eV);
    photonEnergy1.push_back(10324.0*eV);
    photonEnergy1.push_back(4390.1*eV);
    photonEnergy1.push_back(7728.0*eV);
  }
  nphotons = photonEnergy1.size()/2 ;

  for (unsigned i=0; i<nphotons; i++) {
    if (probXaOther > G4UniformRand()) {
      photonEnergy = photonEnergy1[2*i];
    } else {
      photonEnergy = photonEnergy1[2*i+1];
    }
    G4DynamicParticle* producedBoson = new G4DynamicParticle(producedBosonDef, dummy);

//Creation of spearate function to draw in sphere
    G4double costheta = 2.0*G4UniformRand() - 1.0;
    G4double sintheta = std::sqrt(1.0 - costheta*costheta);
    G4double phi = 2.0*pi*G4UniformRand();
    G4double photonMomentum = photonEnergy;

    G4double pz=costheta*photonMomentum;
    G4double px=sintheta*std::cos(phi)*photonMomentum;
    G4double py=sintheta*std::sin(phi)*photonMomentum;

    G4ThreeVector photMomentum(px,py,pz);
    producedBoson->SetMomentum(photMomentum);
    products->push_back(producedBoson);
  }
*/
	// *********************** NEW ******************//
	//

	producedBosonDef = G4Gamma::Gamma();
	G4double photonEnergy;
	G4bool CZT_Al_Target = false;
	G4bool CZT_S_Target = false;
	G4bool CZT_C_Target = false;
	G4bool HPGe_Pb_Target = false;
	unsigned nphotons = 0;

	std::vector<G4double> photonEnergy1;
	if (((iniA) == 27) && ((iniZ) == 13)) {               //Aluminium
		/*  photonEnergy1.push_back(10435.1*eV);
		    photonEnergy1.push_back(17840.0*eV);
		    photonEnergy1.push_back(5113.4*eV);
		    photonEnergy1.push_back(9025.0*eV);

		// added 16.06.2021 H. Shi KAl 7-6 from KG eq. calc.
		photonEnergy1.push_back(16058.0*eV);
		photonEnergy1.push_back(26518.4*eV); // added Jan 2023
		photonEnergy1.push_back(7150.9*eV);
		photonEnergy1.push_back(12262.8*eV);
		 */  } else if (((iniA) == 48) && ((iniZ) == 22)) {        //Titanium
			 photonEnergy1.push_back(8312.0*eV);
			 photonEnergy1.push_back(14770.0*eV);
			 photonEnergy1.push_back(6467.0*eV);
			 photonEnergy1.push_back(11590.0*eV);
		 }/*	else if (((iniA) == 12) && ((iniZ) == 6)) {         //Carbon
			 photonEnergy1.push_back(10216.5*eV);
			 photonEnergy1.push_back(15809.0*eV);
			 photonEnergy1.push_back(5544.9*eV);
			 photonEnergy1.push_back(8885.8*eV);
		 }*/		 else if (((iniA) == 16) && ((iniZ) == 8)) {         //Oxygen
			 photonEnergy1.push_back(9968.7*eV);
			 photonEnergy1.push_back(16062.0*eV);
			 photonEnergy1.push_back(6006.8*eV);
			 photonEnergy1.push_back(9958.0*eV);
		 } else if (((iniA) == 14) && ((iniZ) == 7)) {         //Nitrogen
			 photonEnergy1.push_back(13995.4*eV);
			 photonEnergy1.push_back(21588.2*eV);
			 photonEnergy1.push_back(7595.4*eV);
			 photonEnergy1.push_back(12223*eV);
			 photonEnergy1.push_back(4577.1*eV);
			 photonEnergy1.push_back(7578*eV);
		 } else if (((iniA) == 19) && ((iniZ) == 9)) {         //Fluorine
			 photonEnergy1.push_back(7642.7*eV);
			 photonEnergy1.push_back(12599*eV);
		 } else if (((iniA) == 28) && ((iniZ) == 14)) {        //Silicon
			 photonEnergy1.push_back(8300.1*eV);
			 photonEnergy1.push_back(14233.1*eV);
			 photonEnergy1.push_back(5935.0*eV);
			 photonEnergy1.push_back(10324.0*eV);
			 photonEnergy1.push_back(4390.1*eV);
			 photonEnergy1.push_back(7728.0*eV);
		 }
	nphotons = photonEnergy1.size()/2 ;

	for (unsigned i=0; i<nphotons; i++) {
		if (probXaOther > G4UniformRand()) {
			photonEnergy = photonEnergy1[2*i];
		} else {
			photonEnergy = photonEnergy1[2*i+1];
		}
		G4DynamicParticle* producedBoson = new G4DynamicParticle(producedBosonDef, dummy);

		//Creation of spearate function to draw in sphere
		G4double costheta = 2.0*G4UniformRand() - 1.0;
		G4double sintheta = std::sqrt(1.0 - costheta*costheta);
		G4double phi = 2.0*pi*G4UniformRand();
		G4double photonMomentum = photonEnergy;

		G4double pz=costheta*photonMomentum;
		G4double px=sintheta*std::cos(phi)*photonMomentum;
		G4double py=sintheta*std::sin(phi)*photonMomentum;

		G4ThreeVector photMomentum(px,py,pz);
		producedBoson->SetMomentum(photMomentum);
		products->push_back(producedBoson);
	}
	if (((iniA) == 27) && ((iniZ) == 13))               //Aluminium for CZT
	{
		//G4cout << " IN TARGET ALUMINUM " << G4endl;
		CZT_Al_Target = true;
		photonEnergy1.push_back(302293.0*eV); // 3-->2
		photonEnergy1.push_back(105803.0*eV); // 4-->3
		photonEnergy1.push_back(48972.0*eV); // 5-->4
		photonEnergy1.push_back(75573.0*eV); // 6-->4
		photonEnergy1.push_back(154774.0*eV); // 5-->3
		photonEnergy1.push_back(26602.0*eV); // 6-->5
		photonEnergy1.push_back(42642.0*eV); // 7-->5
	}
	if(CZT_Al_Target)
	{
		nphotons = photonEnergy1.size();

		for (unsigned i=0; i<nphotons; i++) 
		{
			photonEnergy = photonEnergy1[i];
//			G4cout << nphotons << " generated in Al: " << i << " " << photonEnergy / eV << G4endl;

			G4DynamicParticle* producedBoson = new G4DynamicParticle(producedBosonDef, dummy);

			//Creation of spearate function to draw in sphere
			G4double costheta = 2.0*G4UniformRand() - 1.0;
			G4double sintheta = std::sqrt(1.0 - costheta*costheta);
			G4double phi = 2.0*pi*G4UniformRand();
			G4double photonMomentum = photonEnergy;

			G4double pz=costheta*photonMomentum;
			G4double px=sintheta*std::cos(phi)*photonMomentum;
			G4double py=sintheta*std::sin(phi)*photonMomentum;

			G4ThreeVector photMomentum(px,py,pz);
			//	G4cout << px/eV << " " << py/eV << " " << pz/eV << G4endl;
			producedBoson->SetMomentum(photMomentum);
			products->push_back(producedBoson);
		}
		CZT_Al_Target = false;
	}

	if (((iniA) == 12) && ((iniZ) == 6))               //Carbon for CZT
	{
		//G4cout << " IN TARGET Carbon " << G4endl;
		CZT_C_Target = true;
		photonEnergy1.push_back(62881.1*eV); // 3-->2
		photonEnergy1.push_back(22008.4*eV); // 4-->3
		photonEnergy1.push_back(32195.1*eV); // 5-->3
		photonEnergy1.push_back(10216.5*eV); // 5-->4
		photonEnergy1.push_back(15809.0*eV); // 6-->4
		photonEnergy1.push_back(5544.9*eV); // 6-->5
		photonEnergy1.push_back(8885.8*eV); // 7-->5
	}
	if(CZT_C_Target)
	{
		nphotons = photonEnergy1.size();

		for (unsigned i=0; i<nphotons; i++) 
		{
			photonEnergy = photonEnergy1[i];
//			G4cout << nphotons << " generated in C: " << i << " " << photonEnergy / eV << G4endl;

			G4DynamicParticle* producedBoson = new G4DynamicParticle(producedBosonDef, dummy);

			//Creation of spearate function to draw in sphere
			G4double costheta = 2.0*G4UniformRand() - 1.0;
			G4double sintheta = std::sqrt(1.0 - costheta*costheta);
			G4double phi = 2.0*pi*G4UniformRand();
			G4double photonMomentum = photonEnergy;

			G4double pz=costheta*photonMomentum;
			G4double px=sintheta*std::cos(phi)*photonMomentum;
			G4double py=sintheta*std::sin(phi)*photonMomentum;

			G4ThreeVector photMomentum(px,py,pz);
			//	G4cout << px/eV << " " << py/eV << " " << pz/eV << G4endl;
			producedBoson->SetMomentum(photMomentum);
			products->push_back(producedBoson);
		}
		CZT_C_Target = false;
	}

	if (((iniA) == 32) && ((iniZ) == 16))               //Sulfur for CZT
	{
		//G4cout << " IN TARGET SULFUR " << G4endl;
		CZT_S_Target = true;
		photonEnergy1.push_back(160753.0*eV); // 4-->3
		photonEnergy1.push_back(74405.5*eV); // 5-->4
		photonEnergy1.push_back(235158.0*eV); // 5-->3
		photonEnergy1.push_back(40417.8*eV); // 6-->5
		photonEnergy1.push_back(114823.0*eV); // 6-->4
		photonEnergy1.push_back(64788.5*eV); // 7-->5
	}
	if(CZT_S_Target)
	{
		nphotons = photonEnergy1.size();

		for (unsigned i=0; i<nphotons; i++) 
		{
			photonEnergy = photonEnergy1[i];
//			G4cout << nphotons << " generated in Al: " << i << " " << photonEnergy / eV << G4endl;

			G4DynamicParticle* producedBoson = new G4DynamicParticle(producedBosonDef, dummy);

			//Creation of spearate function to draw in sphere
			G4double costheta = 2.0*G4UniformRand() - 1.0;
			G4double sintheta = std::sqrt(1.0 - costheta*costheta);
			G4double phi = 2.0*pi*G4UniformRand();
			G4double photonMomentum = photonEnergy;

			G4double pz=costheta*photonMomentum;
			G4double px=sintheta*std::cos(phi)*photonMomentum;
			G4double py=sintheta*std::sin(phi)*photonMomentum;

			G4ThreeVector photMomentum(px,py,pz);
			//	G4cout << px/eV << " " << py/eV << " " << pz/eV << G4endl;
			producedBoson->SetMomentum(photMomentum);
			products->push_back(producedBoson);
		}
		CZT_S_Target = false;
	}




	if ((((iniA) == 206)||((iniA) == 207)||((iniA) == 208)) && ((iniZ) == 82))               //Lead for CZT
	{
		//G4cout << " IN TARGET ALUMINUM " << G4endl;
		HPGe_Pb_Target = true;
		photonEnergy1.push_back(288822.0*eV); // 9-->8
		photonEnergy1.push_back(206593.0*eV); // 10-->9
		photonEnergy1.push_back(152855.0*eV); // 11-->10
		photonEnergy1.push_back(269114.0*eV); // 12-->10
		photonEnergy1.push_back(359448.0*eV); // 11-->9
		photonEnergy1.push_back(116259.0*eV); // 12-->11
		photonEnergy1.push_back(206736.0*eV); // 13-->11
	}
	if(HPGe_Pb_Target)
	{
		nphotons = photonEnergy1.size();

		for (unsigned i=0; i<nphotons; i++) 
		{
			photonEnergy = photonEnergy1[i];
			//	G4cout << nphotons << " " << i << " " << photonEnergy / eV << G4endl;

			G4DynamicParticle* producedBoson = new G4DynamicParticle(producedBosonDef, dummy);

			//Creation of spearate function to draw in sphere
			G4double costheta = 2.0*G4UniformRand() - 1.0;
			G4double sintheta = std::sqrt(1.0 - costheta*costheta);
			G4double phi = 2.0*pi*G4UniformRand();
			G4double photonMomentum = photonEnergy;

			G4double pz=costheta*photonMomentum;
			G4double px=sintheta*std::cos(phi)*photonMomentum;
			G4double py=sintheta*std::sin(phi)*photonMomentum;

			G4ThreeVector photMomentum(px,py,pz);
			//	G4cout << px/eV << " " << py/eV << " " << pz/eV << G4endl;
			producedBoson->SetMomentum(photMomentum);
			products->push_back(producedBoson);
		}
		HPGe_Pb_Target = false;
	}

//"Special guests": H & D
//generates kH Xray photons only in gas, lot oh H in kapton, Stark killed

  if ((currentMaterial->GetNumberOfElements() == 1) && ((iniA) == 1) && (iniZ) == 1) {
/*
Repeated code -> Create a function and call it!!!
 */
    prob = 0;
    ranflat = G4UniformRand();
    G4double yieldKTot = (yieldKa + yieldKb + yieldKg + yieldKd + yieldKe + yieldKz + yieldKi)/7;

    std::vector<double> photonsEnergyToSim;
    photonsEnergyToSim.push_back(6191.8*eV);
    photonsEnergyToSim.push_back(7388.6*eV);

   /* if ((prob += yieldKa*probKa) > ranflat*yieldKTot) {
      photonEnergy = 6191.8*eV;
    } else if ((prob += yieldKb*probKb) > ranflat*yieldKTot) {
      photonEnergy = 7388.6*eV;
    } else if ((prob += yieldKg*probKg) > ranflat*yieldKTot) {
      photonEnergy = 7807.4*eV;
    } else if ((prob += yieldKd*probKd) > ranflat*yieldKTot) {
      photonEnergy = 8001.3*eV;
    } else if ((prob += yieldKe*probKe) > ranflat*yieldKTot) {
      photonEnergy = 8106.6*eV;
    } else if ((prob += yieldKz*probKz) > ranflat*yieldKTot) {
      photonEnergy = 8170.0*eV;
    } else if ((prob += yieldKi*probKi) > ranflat*yieldKTot) {
      photonEnergy = 8211.2*eV;
    } else {
      G4cout << "======bad probability sum for kH transitions====" << G4endl;
      photonEnergy = 5*eV;
    }*/
    for (unsigned j=0; j<photonsEnergyToSim.size(); j++) {
      photonEnergy = photonsEnergyToSim.at(j);

      if (photonEnergy > 5*eV) {
//uncomment the following line for realistic production (comment the next one)
//G4double photonEnergyH = G4RandBreitWigner::shoot(photonEnergy, 570)*eV;
        G4double photonEnergyH = photonEnergy;
        G4double photonMomentumH = photonEnergyH;

        G4double costhetaH = 2.0*G4UniformRand() - 1.0;
        G4double sinthetaH = std::sqrt(1.0 - costhetaH*costhetaH);
        G4double phiH = 2.0*pi*G4UniformRand();

        G4double pzH = costhetaH*photonMomentumH;
        G4double pxH = sinthetaH*std::cos(phiH)*photonMomentumH;
        G4double pyH = sinthetaH*std::sin(phiH)*photonMomentumH;

        G4ThreeVector photMomentumH(pxH,pyH,pzH);
        G4DynamicParticle* producedBosonH = new G4DynamicParticle(producedBosonDef, dummy);
        producedBosonH->SetMomentum(photMomentumH);
        products->push_back(producedBosonH);
      }
    }
  }

// DEUTERIUM //
//
  if ((currentMaterial->GetNumberOfElements() == 1) && ((iniA) == 2) && (iniZ) == 1) {
    prob = 0;
    ranflat = G4UniformRand();
    G4double shift = -800.*eV;
    G4double yieldKTot = (yieldKa + yieldKb + yieldKg + yieldKd + yieldKe + yieldKz + yieldKi)/7;

    std::vector<double> photonsEnergyToSim;
    photonsEnergyToSim.push_back(7834*eV + shift);
    photonsEnergyToSim.push_back(9280.2*eV + shift);

 /*   if ((prob += yieldKa*probKa) > ranflat*yieldKTot) {
      photonEnergy = 7834.0*eV + shift;
    } else if ((prob += yieldKb*probKb) > ranflat*yieldKTot) {
      photonEnergy = 9280.2*eV + shift;
    } else if ((prob += yieldKg*probKg) > ranflat*yieldKTot) {
      photonEnergy = 9786.2*eV + shift;
    } else if ((prob += yieldKd*probKd) > ranflat*yieldKTot) {
      photonEnergy = 10020.4*eV + shift;
    } else if ((prob += yieldKe*probKe) > ranflat*yieldKTot) {
      photonEnergy = 10147.6*eV + shift;
    } else if ((prob += yieldKz*probKz) > ranflat*yieldKTot) {
      photonEnergy = 10224.3*eV + shift;
    } else if ((prob += yieldKi*probKi) > ranflat*yieldKTot) {
      photonEnergy = 10274.1*eV + shift;
    } else {
      G4cout << "======bad probability sum for kH transitions====" << G4endl;
      photonEnergy = 5*eV;
    }*/
    for (unsigned j=0; j<photonsEnergyToSim.size(); j++) {
      photonEnergy = photonsEnergyToSim.at(j);

      if (photonEnergy > 5*eV) {

//uncomment the following line for realistic production (comment the next one)
//G4double photonEnergyH = G4RandBreitWigner::shoot(photonEnergy, 570)*eV;
        G4double photonEnergyH = photonEnergy;
        G4double photonMomentumH = photonEnergyH;

        G4double costhetaH = 2.0*G4UniformRand() - 1.0;
        G4double sinthetaH = std::sqrt(1.0 - costhetaH*costhetaH);
        G4double phiH = 2.0*pi*G4UniformRand();

        G4double pzH=costhetaH*photonMomentumH;
        G4double pxH=sinthetaH*std::cos(phiH)*photonMomentumH;
        G4double pyH=sinthetaH*std::sin(phiH)*photonMomentumH;

        G4ThreeVector photMomentumH(pxH,pyH,pzH);
        G4DynamicParticle* producedBosonH = new G4DynamicParticle(producedBosonDef,dummy);
        producedBosonH->SetMomentum(photMomentumH);
        products->push_back(producedBosonH);
      }
    }
  }

// Helium-4 // Added 12.2020 H. Shi
//
  if ((currentMaterial->GetNumberOfElements() == 1) && ((iniA) == 4) && (iniZ) == 2) {
// Assume 100 % yield La line only
// shift from PLB 681 310
    G4double shift = 0.*eV;
    photonEnergy = 6463.6*eV + shift;

    if (photonEnergy > 5*eV) {
//uncomment the following line for realistic production (comment the next one)
//G4double photonEnergyH = G4RandBreitWigner::shoot(photonEnergy, 570)*eV;
      G4double photonEnergyH = photonEnergy;
      G4double photonMomentumH = photonEnergyH;

      G4double costhetaH = 2.0*G4UniformRand() - 1.0;
      G4double sinthetaH = std::sqrt(1.0 - costhetaH*costhetaH);
      G4double phiH = 2.0*pi*G4UniformRand();

      G4double pzH=costhetaH*photonMomentumH;
      G4double pxH=sinthetaH*std::cos(phiH)*photonMomentumH;
      G4double pyH=sinthetaH*std::sin(phiH)*photonMomentumH;

      G4ThreeVector photMomentumH(pxH,pyH,pzH);
      G4DynamicParticle* producedBosonH = new G4DynamicParticle(producedBosonDef, dummy);
      producedBosonH->SetMomentum(photMomentumH);
      products->push_back(producedBosonH);
    }
  }

//Neonium
  if ((currentMaterial->GetNumberOfElements() == 1) && ((iniA) == 20) && (iniZ) == 10) {
// Assume 100 % yield La line only
    G4double shift = 0.*eV;
//Parameter to choose the line to simulate, instead of saving commented lines. In future, maybe it is worth to manipulate the number
//of line by modifying CARD.dat or using messenger.
    int whichLine = 1; // 2 -> second , 3 -> third
    std::vector<G4double> photonEnergyLines;

    photonEnergyLines.push_back(6118.91*eV + shift); //8->7
    photonEnergyLines.push_back(9427.65*eV + shift); //7->6
    photonEnergyLines.push_back(15635.4*eV + shift); //6->5

// Simulating only one line at a time
//photonEnergy = photonEnergyLines.at(whichLine - 1);

// Simulating all three lines at once
    for (unsigned j=0; j<photonEnergyLines.size(); j++) {

      photonEnergy = photonEnergyLines.at(j);

      if (photonEnergy > 5*eV) {
//uncomment the following line for realistic production (comment the next one)
//G4double photonEnergyH = G4RandBreitWigner::shoot(photonEnergy, 570)*eV;
        G4double photonEnergyH = photonEnergy; 
       	G4double photonMomentumH = photonEnergyH;

	G4double costhetaH = 2.0*G4UniformRand() - 1.0;
	G4double sinthetaH = std::sqrt(1.0 - costhetaH*costhetaH);
	G4double phiH = 2.0*pi*G4UniformRand();

	G4double pzH=costhetaH*photonMomentumH;
	G4double pxH=sinthetaH*std::cos(phiH)*photonMomentumH;
	G4double pyH=sinthetaH*std::sin(phiH)*photonMomentumH;

	G4ThreeVector photMomentumH(pxH,pyH,pzH);
	G4DynamicParticle* producedBosonH = new
	G4DynamicParticle(producedBosonDef, dummy);
	producedBosonH->SetMomentum(photMomentumH);
	products->push_back(producedBosonH);
      }
    }
  } 
//==========================end cascade==========================================
  if (fPlotProducedXrayLies) {
    SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
    for (unsigned it=0; it<products->size(); it++) {
      analysis->histo->fillHistogram("kaonAbsorptionEnergies_vs_Nucleus", products->at(it)->GetTotalMomentum()/eV, doubleCheck(iniA), doubleCheck(iniA));
    }
  }

  products->push_back(producedBaryon);
  products->push_back(producedMeson);

  if (verboseLevel > 1) {
    G4cout << "G4KaonMinusAbsorption::KaonNucleonReaction: Number of primaries = " << products->size()
           << ": " << producedMesonDef->GetParticleName() << ", " << producedBaryonDef->GetParticleName()
           << ", " << producedBosonDef->GetParticleName() << G4endl;
  }
  return products;
}


G4bool G4KaonMinusAbsorptionAtRest::AbsorbPionByNucleus(G4DynamicParticle* aPion)
{
// Needs some more investigation!
  G4double ranflat = G4UniformRand();

  if (ranflat < pionAbsorptionRate) {
// Add pion energy to ExcitationEnergy and NucleusMomentum
    nucleus->AddExcitationEnergy(aPion->GetTotalEnergy());
    nucleus->AddMomentum(aPion->GetMomentum());
  }

  return (ranflat < pionAbsorptionRate);
}


G4DynamicParticle* G4KaonMinusAbsorptionAtRest::SigmaLambdaConversion(G4DynamicParticle* aSigma)
{
  G4double  ranflat = G4UniformRand();
  G4double  sigmaLambdaConversionRate;

  G4int A = nucleus->GetA_asInt();
  G4int Z = nucleus->GetZ_asInt();

  G4int newZ = Z;
  G4double nucleonMassDifference = 0;

  G4ParticleDefinition* inNucleonDef = nullptr;
  G4ParticleDefinition* outNucleonDef = nullptr;

  // Decide which sigma
  switch((int)aSigma->GetDefinition()->GetPDGCharge()) {
    case 1:
      sigmaLambdaConversionRate = sigmaPlusLambdaConversionRate;
      inNucleonDef   = G4Neutron::Neutron();
      outNucleonDef  = G4Proton::Proton();
      newZ = Z+1;
      nucleonMassDifference = neutron_mass_c2 - (proton_mass_c2 + electron_mass_c2);
      break;
    case -1:
      sigmaLambdaConversionRate = sigmaMinusLambdaConversionRate;
      inNucleonDef   = G4Proton::Proton();
      outNucleonDef  = G4Neutron::Neutron();
      newZ = Z-1;
      nucleonMassDifference = (proton_mass_c2 + electron_mass_c2) - neutron_mass_c2;
      break;
    case 0:
      sigmaLambdaConversionRate = sigmaZeroLambdaConversionRate;
// The 'outgoing' nucleon is just virtual, to keep the energy-momentum
// balance and will not appear in the ParticleChange. Therefore no need
// choose between neutron and proton here!
      inNucleonDef   = G4Neutron::Neutron();
      outNucleonDef  = G4Neutron::Neutron();
      break;
    default:
      sigmaLambdaConversionRate = 0.;
  }
  if (ranflat >= sigmaLambdaConversionRate)
    return nullptr;

  G4ThreeVector dummy(0.,0.,0.);

// Fermi momentum distribution in three dimensions
  G4ThreeVector momentum = nucleus->GetFermiMomentum();
  G4ParticleDefinition* lambdaDef  = G4Lambda::Lambda();

  G4DynamicParticle inNucleon(inNucleonDef, momentum);
  G4DynamicParticle outNucleon(outNucleonDef, dummy);
  G4DynamicParticle* outLambda = new G4DynamicParticle(lambdaDef, dummy);
  G4ReactionKinematics theReactionKinematics;
// Now do the twobody scattering
  theReactionKinematics.TwoBodyScattering(aSigma, &inNucleon, &outNucleon, outLambda);

// Binding energy of nucleus has changed. This will change the
// ExcitationEnergy.
// .. using mass formula ..
// .. using mass table ..
// equivalent to -'initialBindingEnergy+nucleus.GetBindingEnergy' !

// Add energy and momentum to nucleus, change Z,A
  nucleus->AddExcitationEnergy(outNucleon.GetKineticEnergy());
  nucleus->AddMomentum(outNucleon.GetMomentum());

  if (newZ<0) //Problem for deuterium -> what is a neutron electron pair?!
    newZ = 0;
  nucleus->SetParameters(A,newZ);
// The calling routine is responsible to delete the sigma!!
  return outLambda;
}
