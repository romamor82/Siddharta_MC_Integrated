#include "../include/SiddhartaPhysicsList.h"
#include "../include/SiddhartaCard.h"

#include "../include/SpecialCuts.h"

#include <G4ShortLivedConstructor.hh>
#include <G4ParticleDefinition.hh>
#include <G4LeptonConstructor.hh>
#include <G4BaryonConstructor.hh>
#include <G4MesonConstructor.hh>
#include <G4ProcessManager.hh>
#include <G4IonConstructor.hh>
#include <G4ParticleTypes.hh>
#include <G4EmParameters.hh>
#include <G4DecayWithSpin.hh>
#include <globals.hh>

#include <G4ParticleWithCuts.hh>
#include <G4ParticleTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4ProcessVector.hh>
#include <G4UserLimits.hh>
#include <G4ios.hh>

SiddhartaPhysicsList::SiddhartaPhysicsList():  G4VUserPhysicsList()
{
  SiddhartaCard* mycard = SiddhartaCard::getInstance();
  int SiddhartaSetup = mycard->variables["SiddhartaSetupVersion"];
  defaultCutValue = mycard->variables["defaultCutValue"]*micrometer;
  LowEnCut = mycard->variables["lowlimit"]*eV;
  cutForGamma 	  = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
  cutForKaonPlus  = defaultCutValue;
  cutForKaonMinus = defaultCutValue;
  cutForPionPlus  = defaultCutValue;
  cutForPionMinus = defaultCutValue;
  cutForProton    = defaultCutValue;
  cutForNeutron   = defaultCutValue;
  SetVerboseLevel(0);

  G4EmParameters* param = G4EmParameters::Instance();
  param->SetMaxEnergy(100*GeV);
  param->SetNumberOfBinsPerDecade(20);
  param->SetMscStepLimitType(fMinimal);
  param->SetFluo(true);
  param->SetPixe(true);
  param->SetAuger(true);
  param->SetDeexActiveRegion("World" , true, true, true);
}

SiddhartaPhysicsList::~SiddhartaPhysicsList() {}

void SiddhartaPhysicsList::ConstructBosons()
{
//pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  G4Gamma::GammaDefinition();
}

void SiddhartaPhysicsList::ConstructLeptons()
{
  G4LeptonConstructor lConstructor;
  lConstructor.ConstructParticle();
}

void SiddhartaPhysicsList::ConstructMesons()
{
  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();
}


 void SiddhartaPhysicsList::ConstructShortLived()
{
  G4ShortLivedConstructor shor;
  shor.ConstructParticle();
}

void SiddhartaPhysicsList::ConstructBaryons()
{
  G4BaryonConstructor baryon;
  baryon.ConstructParticle();
}


void SiddhartaPhysicsList::ConstructIons()
{
  G4IonConstructor ions;
  ions.ConstructParticle();
}

#include "G4StepLimiter.hh"

void SiddhartaPhysicsList::AddTransportation()
{
  G4VUserPhysicsList::AddTransportation();
  auto particleIterator = GetParticleIterator();
  particleIterator->reset();

  while ((*particleIterator)()) {
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if(particleName == "neutron")
      pmanager->AddDiscreteProcess(new MaxTimeCuts());
    pmanager->AddDiscreteProcess(new MinEkineCuts());

    pmanager->AddProcess(new G4StepLimiter, -1, -1, 1);
  }
}

#include <G4PhotoElectricEffect.hh>
#include <G4BetheHeitler5DModel.hh>
#include <G4ComptonScattering.hh>
#include <G4GammaConversion.hh>

#include <G4eMultipleScattering.hh>
#include <G4eplusAnnihilation.hh>
#include <G4eBremsstrahlung.hh>
#include <G4eIonisation.hh>

#include <G4MuMultipleScattering.hh>
#include <G4MuPairProduction.hh>
#include <G4MuBremsstrahlung.hh>
#include <G4MuonMinusCapture.hh>
#include <G4MuIonisation.hh>

#include <G4hMultipleScattering.hh>
#include <G4hPairProduction.hh>
#include <G4hBremsstrahlung.hh>
#include <G4hIonisation.hh>

#include <G4IonParametrisedLossModel.hh>
#include <G4ionIonisation.hh>

#include <G4LivermoreGammaConversionModel.hh>
#include <G4LivermoreBremsstrahlungModel.hh>
#include <G4LivermorePhotoElectricModel.hh>
#include <G4LivermoreIonisationModel.hh>
#include <G4LivermoreRayleighModel.hh>
#include <G4LivermoreComptonModel.hh>
#include <G4UniversalFluctuation.hh>
#include <G4RayleighScattering.hh>
#include <G4LossTableManager.hh>
#include <G4EnergyLossTables.hh>
//#include <G4EmProcessOptions.hh>

#include <G4UAtomicDeexcitation.hh>
#include <G4VAtomDeexcitation.hh>

void SiddhartaPhysicsList::ConstructEM()
{
  G4LossTableManager* man = G4LossTableManager::Instance();
  G4VAtomDeexcitation* ad = man->AtomDeexcitation();

  if(!ad) {
    man->SetAtomDeexcitation(new G4UAtomicDeexcitation());
  }

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();

  while((*particleIterator)()){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    G4String particleType = particle->GetParticleType();
    G4double charge = particle->GetPDGCharge();

    if (particleName == "gamma") {
      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      pmanager->AddDiscreteProcess(theRayleigh);

      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel());
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);

      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      theComptonScattering->SetEmModel(new G4LivermoreComptonModel());
      pmanager->AddDiscreteProcess(theComptonScattering);

      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      theGammaConversion->SetEmModel(new G4BetheHeitler5DModel());
      pmanager->AddDiscreteProcess(theGammaConversion);
    } else if (particleName == "e-") {
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc, -1, 1, 1);

      G4eIonisation* eIonisation = new G4eIonisation();
      G4VEmModel* theIoniLiv = new G4LivermoreIonisationModel();
      theIoniLiv->SetHighEnergyLimit(0.1*MeV);
      eIonisation->AddEmModel(0, theIoniLiv, new G4UniversalFluctuation() );
      eIonisation->SetStepFunction(0.2, 100*um); //improved precision in tracking
      pmanager->AddProcess(eIonisation, -1, 2, 2);

      G4eBremsstrahlung* eBremsstrahlung = new G4eBremsstrahlung();
      pmanager->AddProcess(eBremsstrahlung, -1, -3, 3);
    } else if (particleName == "e+") {
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc, -1, 1, 1);

      G4eIonisation* eIonisation = new G4eIonisation();
      eIonisation->SetStepFunction(0.2, 100*um); //
      pmanager->AddProcess(eIonisation, -1, 2, 2);

      pmanager->AddProcess(new G4eBremsstrahlung(), -1, -3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation(), 0, -1, 4);
    } else if(particleName == "mu+" || particleName == "mu-") {
      pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation(), -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung(), -1, 3, 3);
      pmanager->AddProcess(new G4MuPairProduction(), -1, 4, 4);
      if(particleName == "mu-")
        pmanager->AddProcess(new G4MuonMinusCapture(), 0, -1, -1);
    } else if (particleName == "proton" || particleName == "pi-" || particleName == "pi+"  || particleName == "kaon-" || particleName == "kaon+") {
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      G4hIonisation* hIonisation = new G4hIonisation();
     // hIonisation->SetStepFunction(0.2, 50*um);
      pmanager->AddProcess(hIonisation, -1, 2, 2);
      pmanager->AddProcess(new G4hBremsstrahlung, -1, 3, 3);
      pmanager->AddProcess(new G4hPairProduction, -1, 4, 4);
    } else if(particleName == "alpha" || particleName == "He3" || particleName == "GenericIon") {
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetStepFunction(0.1, 20*um);
      pmanager->AddProcess(ionIoni, -1, 2, 2);
    } else if ((!particle->IsShortLived()) && (charge != 0.0) && (particle->GetParticleName() != "chargedgeantino")) {
      G4hMultipleScattering* aMultipleScattering = new G4hMultipleScattering();
      G4hIonisation* ahadronIon = new G4hIonisation();
	  pmanager->AddProcess(aMultipleScattering, -1, 1, 1);
	  pmanager->AddProcess(ahadronIon, -1, 2, 2);
    }
  }
}

//========================Hadronic processes==============================

// Elastic processes:
#include <G4ElasticHadrNucleusHE.hh>
#include <G4HadronElasticProcess.hh>
#include <G4ChipsElasticModel.hh>
#include <G4HadronElastic.hh>

// Inelastic processes:
#include <G4HadronInelasticProcess.hh>

// High energy FTFP model and Bertini cascade
#include <G4GeneratorPrecompoundInterface.hh>
#include <G4LundStringFragmentation.hh>
#include <G4ExcitedStringDecay.hh>
#include <G4PreCompoundModel.hh>
#include <G4CascadeInterface.hh>
#include <G4TheoFSGenerator.hh>
#include <G4FTFModel.hh>

// Cross sections
#include <G4CrossSectionDataSetRegistry.hh>
#include <G4VCrossSectionDataSet.hh>

#include <G4CrossSectionInelastic.hh>
#include <G4CrossSectionElastic.hh>
#include <G4BGGPionInelasticXS.hh>
#include <G4BGGPionElasticXS.hh>
#include <G4AntiNuclElastic.hh>

#include <G4ComponentGGHadronNucleusXsc.hh>
#include <G4ComponentAntiNuclNuclearXS.hh>
#include <G4ComponentGGNuclNuclXsc.hh>
#include <G4CrossSectionInelastic.hh>
#include <G4BGGNucleonInelasticXS.hh>
#include <G4BGGNucleonElasticXS.hh>
#include <G4NeutronInelasticXS.hh>
#include <G4NeutronElasticXS.hh>

#include <G4NeutronCaptureProcess.hh>
#include <G4HadronElastic.hh>

// Neutron high-precision models: <20 MeV
#include <G4ParticleHPInelasticData.hh>
#include <G4ParticleHPElasticData.hh>
#include <G4ParticleHPCaptureData.hh>
#include <G4ParticleHPInelastic.hh>
#include <G4ParticleHPCapture.hh>
#include <G4ParticleHPElastic.hh>

// Stopping processes
#include "../include/G4KaonMinusAbsorptionAtRest.h"

//#include <G4AntiProtonAbsorptionFritiof.hh>
//#include <G4KaonMinusAbsorptionBertini.hh>
//#include <G4PiMinusAbsorptionBertini.hh>
#include <G4HadronicParameters.hh>
#include <G4HadronicAbsorptionBertini.hh>

void SiddhartaPhysicsList::ConstructHadronic()
{
  const G4double elastic_elimitPi = 1.0*GeV;

  G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE();

  // Inelastic scattering
  const G4double theFTFMin0 = 0.0*GeV;
  const G4double theFTFMin1 = 4.0*GeV;
  const G4double theFTFMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  const G4double theBERTMin0 = 0.0*GeV;
  const G4double theBERTMin1 = 19.0*MeV;
  const G4double theBERTMax = 6.0*GeV;
  const G4double theHPMin = 0.0*GeV;
  const G4double theHPMax = 20.0*MeV;

  G4FTFModel * theStringModel = new G4FTFModel;
  G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay(new G4LundStringFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);
  G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
  G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface(thePreEquilib);

  G4TheoFSGenerator * theFTFModel0 = new G4TheoFSGenerator("FTFP");
  theFTFModel0->SetHighEnergyGenerator(theStringModel);
  theFTFModel0->SetTransport(theCascade);
  theFTFModel0->SetMinEnergy(theFTFMin0);
  theFTFModel0->SetMaxEnergy(theFTFMax);

  G4TheoFSGenerator * theFTFModel1 = new G4TheoFSGenerator("FTFP");
  theFTFModel1->SetHighEnergyGenerator(theStringModel);
  theFTFModel1->SetTransport(theCascade);
  theFTFModel1->SetMinEnergy(theFTFMin1);
  theFTFModel1->SetMaxEnergy(theFTFMax);

  G4CascadeInterface * theBERTModel0 = new G4CascadeInterface;
  theBERTModel0->SetMinEnergy(theBERTMin0);
  theBERTModel0->SetMaxEnergy(theBERTMax);

  G4CascadeInterface * theBERTModel1 = new G4CascadeInterface;
  theBERTModel1->SetMinEnergy(theBERTMin1);
  theBERTModel1->SetMaxEnergy(theBERTMax);

  G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic(new G4ComponentAntiNuclNuclearXS);
  G4ComponentGGNuclNuclXsc * ggNuclNuclXsec = new G4ComponentGGNuclNuclXsc();
  G4VCrossSectionDataSet * theGGNuclNuclData = new G4CrossSectionInelastic(ggNuclNuclXsec);
  G4VCrossSectionDataSet * theGGNNEl = new G4CrossSectionElastic(ggNuclNuclXsec);
  G4ComponentGGHadronNucleusXsc * ggHNXsec = new G4ComponentGGHadronNucleusXsc();
  G4VCrossSectionDataSet * theGGHNEl = new G4CrossSectionElastic(ggHNXsec);
  G4VCrossSectionDataSet * theGGHNInel = new G4CrossSectionInelastic(ggHNXsec);

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)()) {
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "pi+") {
// Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(new G4BGGPionElasticXS(particle));
      theElasticProcess->RegisterMe(elastic_he);
      pmanager->AddDiscreteProcess(theElasticProcess);
//Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(new G4BGGPionInelasticXS(particle));
	  theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    } else if (particleName == "pi-") {
// Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(new G4BGGPionElasticXS(particle));
      theElasticProcess->RegisterMe(elastic_he);
      pmanager->AddDiscreteProcess(theElasticProcess);
//Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(new G4BGGPionInelasticXS(particle));
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
//Absorption
  //    pmanager->AddRestProcess(new G4PiMinusAbsorptionBertini, ordDefault);
    } else if (particleName == "kaon+") {
// Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGHNEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
// Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(theGGHNInel);
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    } else if (particleName == "kaon0S") {
// Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGHNEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
// Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(theGGHNInel);
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    } else if (particleName == "kaon0L") {
// Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGHNEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
// Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(theGGHNInel);
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    } else if (particleName == "kaon-") {
// Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGHNEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
// Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(theGGHNInel);
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4KaonMinusAbsorptionAtRest, ordDefault);    
    } else if (particleName == "proton") {
// Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(new G4BGGNucleonElasticXS(G4Proton::Proton()));
      theElasticProcess->RegisterMe(elastic_chip);
      pmanager->AddDiscreteProcess(theElasticProcess);
// Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(new G4BGGNucleonInelasticXS(G4Proton::Proton()));
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    } else if (particleName == "anti_proton") {
// Elastic scattering
      const G4double elastic_elimitAntiNuc = 100.0*MeV;
      G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
      elastic_anuc->SetMinEnergy(elastic_elimitAntiNuc);
      G4CrossSectionElastic* elastic_anucxs = new G4CrossSectionElastic(elastic_anuc->GetComponentCrossSection());
      G4HadronElastic* elastic_lhep2 = new G4HadronElastic();
      elastic_lhep2->SetMaxEnergy(elastic_elimitAntiNuc);
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(elastic_anucxs);
      theElasticProcess->RegisterMe(elastic_lhep2);
      theElasticProcess->RegisterMe(elastic_anuc);
      pmanager->AddDiscreteProcess(theElasticProcess);
// Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(theAntiNucleonData);
      theInelasticProcess->RegisterMe(theFTFModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
// Absorption
      pmanager->AddRestProcess(new G4HadronStoppingProcess, ordDefault);
    } else if (particleName == "neutron") {
// elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(new G4NeutronElasticXS());
      G4HadronElastic* elastic_neutronChipsModel = new G4ChipsElasticModel();
      elastic_neutronChipsModel->SetMinEnergy(19.0*MeV);
      theElasticProcess->RegisterMe(elastic_neutronChipsModel);
      G4ParticleHPElastic * theElasticNeutronHP = new G4ParticleHPElastic;
      theElasticNeutronHP->SetMinEnergy(theHPMin);
      theElasticNeutronHP->SetMaxEnergy(theHPMax);
      theElasticProcess->RegisterMe(theElasticNeutronHP);
      theElasticProcess->AddDataSet(new G4ParticleHPElasticData);
      pmanager->AddDiscreteProcess(theElasticProcess);
// inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(new G4NeutronInelasticXS());
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel1);
      G4ParticleHPInelastic * theNeutronInelasticHPModel = new G4ParticleHPInelastic;
      theNeutronInelasticHPModel->SetMinEnergy(theHPMin);
      theNeutronInelasticHPModel->SetMaxEnergy(theHPMax);
      theInelasticProcess->RegisterMe(theNeutronInelasticHPModel);
      theInelasticProcess->AddDataSet(new G4ParticleHPInelasticData);
      pmanager->AddDiscreteProcess(theInelasticProcess);
// capture
      G4NeutronCaptureProcess* theCaptureProcess = new G4NeutronCaptureProcess;
      G4ParticleHPCapture * theLENeutronCaptureModel = new G4ParticleHPCapture;
      theLENeutronCaptureModel->SetMinEnergy(theHPMin);
      theLENeutronCaptureModel->SetMaxEnergy(theHPMax);
      theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
      theCaptureProcess->AddDataSet(new G4ParticleHPCaptureData);
      pmanager->AddDiscreteProcess(theCaptureProcess);
    } else if (particleName == "anti_neutron") {
// Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGHNEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
// Inelastic scattering (include annihilation on-fly)
      G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(theAntiNucleonData);
      theInelasticProcess->RegisterMe(theFTFModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    } else if (particleName == "deuteron") {
// Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGNNEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
// Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(theGGNuclNuclData);
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    } else if (particleName == "triton") {
// Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGNNEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
// Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(theGGNuclNuclData);
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    } else if (particleName == "alpha") {
// Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGNNEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
// Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(theGGNuclNuclData);
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
}

//==================================================================
#include <G4PhysicsListHelper.hh>
#include <G4RadioactiveDecay.hh>
#include <G4Decay.hh>

void SiddhartaPhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  auto theParticleIterator = GetParticleIterator();
  theParticleIterator->reset();

  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (theDecayProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(theDecayProcess);
      pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);

      if (particleName != "kaon-")
        pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }

  G4LossTableManager* man = G4LossTableManager::Instance();
  G4VAtomDeexcitation* ad = man->AtomDeexcitation();
  if(!ad) {
    G4EmParameters::Instance()->SetAugerCascade(true);
    ad = new G4UAtomicDeexcitation();
    man->SetAtomDeexcitation(ad);
    ad->InitialiseAtomicDeexcitation();
  }

  G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(new G4RadioactiveDecay(), G4GenericIon::GenericIon());
}

#include <G4UserSpecialCuts.hh>

void SiddhartaPhysicsList::AddStepMax()
{
  G4StepLimiter* stepLimiter = new G4StepLimiter();
  auto theParticleIterator = GetParticleIterator();
  theParticleIterator->reset();

  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    
    if (particle->GetPDGCharge() != 0.0) {
      pmanager->AddDiscreteProcess(stepLimiter);
    }
  }
}

//========================Low Energy Cuts=========================

void SiddhartaPhysicsList::SetCuts()
{
  if (verboseLevel>0) {
    G4cout << "SiddhartaPhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(LowEnCut,100.*GeV);

  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  SetCutValue(cutForKaonPlus, "kaon+");
  SetCutValue(cutForKaonMinus,"kaon-");
  SetCutValue(cutForPionPlus, "pi+");
  SetCutValue(cutForPionMinus,"pi-");
  SetCutValue(cutForProton, "proton");
  SetCutValue(cutForNeutron,"neutron");

  if (verboseLevel>0)
    DumpCutValuesTable();
}

void SiddhartaPhysicsList::ConstructParticle()
{
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructShortLived();
  ConstructIons();
}

void SiddhartaPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructHadronic();
  ConstructGeneral();
  AddStepMax();
  SetCuts();
}
