#include "../include/SinsqDecayChannel.h"

#include <G4ParticleDefinition.hh>
#include <G4PhysicalConstants.hh>
#include <G4LorentzRotation.hh>
#include <G4VDecayChannel.hh>
#include <G4LorentzVector.hh>
#include <G4DecayProducts.hh>
#include <Randomize.hh>

SinsqDecayChannel::SinsqDecayChannel(G4int Verbose) : G4VDecayChannel("Phase Space", Verbose) {}

SinsqDecayChannel::SinsqDecayChannel(const G4String& theParentName, G4double theBR, G4int theNumberOfDaughters, const G4String& theDaughterName1,
                                     const G4String& theDaughterName2, const G4String& theDaughterName3, const G4String& theDaughterName4)
                         : G4VDecayChannel("Phase Space", theParentName, theBR, theNumberOfDaughters,
                                           theDaughterName1, theDaughterName2, theDaughterName3, theDaughterName4) {}

SinsqDecayChannel::~SinsqDecayChannel() {}

G4DecayProducts *SinsqDecayChannel::DecayIt(G4double parentMass)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1)
    G4cout << "SinsqDecayChannel::DecayIt ";
#endif

  G4DecayProducts * products = 0;

  if (G4MT_parent == 0)
    CheckAndFillParent();

  if (G4MT_daughters == 0)
    CheckAndFillDaughters();

  if (G4MT_parent_mass > 0.0)
    G4MT_parent_mass = parentMass;

  switch (numberOfDaughters) {
    case 0:
#ifdef G4VERBOSE
      if (GetVerboseLevel() > 0) {
        G4cout << "SinsqDecayChannel::DecayIt daughters not defined " << G4endl;
      }
#endif
      break;
    case 1:
      products = OneBodyDecayIt();
      break;
    case 2:
      products = TwoBodyDecayIt();
      break;
    case 3:
      products = ThreeBodyDecayIt();
      break;
    default:
      products = ManyBodyDecayIt();
      break;
  }
#ifdef G4VERBOSE
  if ((products == 0) && (GetVerboseLevel()>0)) {
    G4cout << "SinsqDecayChannel::DecayIt " << *parent_name << " can not decay " << G4endl;
    DumpInfo();
  }
#endif
  return products;
}

G4DecayProducts *SinsqDecayChannel::OneBodyDecayIt()
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1)
    G4cout << "SinsqDecayChannel::OneBodyDecayIt()" << G4endl;
#endif

  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle(G4MT_parent, dummy, 0.0);
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  G4DynamicParticle * daughterparticle = new G4DynamicParticle(G4MT_daughters[0], dummy, 0.0);
  products->PushProducts(daughterparticle);

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
     G4cout << "SinsqDecayChannel::OneBodyDecayIt create decay products in rest frame " << G4endl;
     products->DumpInfo();
  }
#endif
  return products;
}

G4DecayProducts *SinsqDecayChannel::TwoBodyDecayIt()
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1)
    G4cout << "SinsqDecayChannel::TwoBodyDecayIt()" << G4endl;
#endif

  G4double parentmass = G4MT_parent_mass;
  G4double daughtermass[2];
  G4double daughtermomentum;
  daughtermass[0] = G4MT_daughters_mass[0];
  daughtermass[1] = G4MT_daughters_mass[1];

  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle(G4MT_parent, dummy, 0.0);
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  daughtermomentum = Pmx(parentmass, daughtermass[0], daughtermass[1]);
  if (daughtermomentum < 0.0) {
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 0) {
      G4cerr << "SinsqDecayChannel::TwoBodyDecayIt " << "sum of daughter mass is larger than parent mass" << G4endl;
      G4cerr << "parent :" << G4MT_parent->GetParticleName() << "  " << G4MT_parent_mass/GeV << G4endl;
      G4cerr << "daughter 1 :" << G4MT_daughters[0]->GetParticleName() << "  " << daughtermass[0]/GeV << G4endl;
      G4cerr << "daughter 2:" << G4MT_daughters[1]->GetParticleName() << "  " << daughtermass[1]/GeV << G4endl;
    }
#endif
    G4Exception("SinsqDecayChannel::TwoBodyDecayIt", "can not create decay products", JustWarning, "sum of daughter mass is larger than parent mass");
    return products;
  }
  G4bool condition = 0;
  G4double costheta;
  G4double sintheta;

  do {
   costheta = 2.*G4UniformRand() - 1.0;
   sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta));
   G4double cut = G4UniformRand();
   condition = sintheta*sintheta >= cut;
  } while (!condition);
  G4double phi = twopi*G4UniformRand()*rad;
  G4ThreeVector direction(sintheta*std::cos(phi), sintheta*std::sin(phi), costheta);

  G4DynamicParticle * daughterparticle = new G4DynamicParticle(G4MT_daughters[0], direction*daughtermomentum);
  products->PushProducts(daughterparticle);
  daughterparticle = new G4DynamicParticle(G4MT_daughters[1], direction*(-1.0*daughtermomentum));
  products->PushProducts(daughterparticle);

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
     G4cout << "SinsqDecayChannel::TwoBodyDecayIt create decay products in rest frame " << G4endl;
     products->DumpInfo();
  }
#endif
  return products;
}

G4DecayProducts *SinsqDecayChannel::ThreeBodyDecayIt()
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1)
    G4cout << "SinsqDecayChannel::ThreeBodyDecayIt()" << G4endl;
#endif

  G4double parentmass = G4MT_parent_mass;
  G4double daughtermass[3];
  G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<3; index++) {
    daughtermass[index] = G4MT_daughters_mass[index];
    sumofdaughtermass += daughtermass[index];
  }
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle(G4MT_parent, dummy, 0.0);
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  if (sumofdaughtermass > parentmass) {
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 0) {
      G4cerr << "SinsqDecayChannel::ThreeBodyDecayIt " << "sum of daughter mass is larger than parent mass" << G4endl;
      G4cerr << "parent :" << G4MT_parent->GetParticleName() << "  " << G4MT_parent_mass/GeV << G4endl;
      G4cerr << "daughter 1 :" << G4MT_daughters[0]->GetParticleName() << "  " << daughtermass[0]/GeV << G4endl;
      G4cerr << "daughter 2:" << G4MT_daughters[1]->GetParticleName() << "  " << daughtermass[1]/GeV << G4endl;
      G4cerr << "daughter 3:" << G4MT_daughters[2]->GetParticleName() << "  " << daughtermass[2]/GeV << G4endl;
    }
#endif
    G4Exception("SinsqDecayChannel::ThreeBodyDecayIt", "can not create decay products", JustWarning, "sum of daughter mass is larger than parent mass");
    return products;
  }
  G4double rd1, rd2, rd;
  G4double daughtermomentum[3];
  G4double momentummax = 0.0, momentumsum = 0.0;
  G4double energy;

  do {
    rd1 = G4UniformRand();
    rd2 = G4UniformRand();

    if (rd2 > rd1) {
      rd = rd1;
      rd1 = rd2;
      rd2 = rd;
    }
    momentummax = 0.0;
    momentumsum = 0.0;
    energy = rd2*(parentmass - sumofdaughtermass);
    daughtermomentum[0] = std::sqrt(energy*energy + 2.0*energy* daughtermass[0]);

    if (daughtermomentum[0] > momentummax)
      momentummax = daughtermomentum[0];
    momentumsum += daughtermomentum[0];
    energy = (1. - rd1)*(parentmass - sumofdaughtermass);
    daughtermomentum[1] = std::sqrt(energy*energy + 2.0*energy*daughtermass[1]);

    if (daughtermomentum[1] > momentummax)
      momentummax = daughtermomentum[1];
    momentumsum += daughtermomentum[1];
    energy = (rd1 - rd2)*(parentmass - sumofdaughtermass);
    daughtermomentum[2] = std::sqrt(energy*energy + 2.0*energy* daughtermass[2]);

    if (daughtermomentum[2] > momentummax)
      momentummax = daughtermomentum[2];
    momentumsum += daughtermomentum[2];
  } while (momentummax > momentumsum - momentummax);
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "     daughter 0:" << daughtermomentum[0]/GeV << "[GeV/c]" << G4endl;
    G4cout << "     daughter 1:" << daughtermomentum[1]/GeV << "[GeV/c]" << G4endl;
    G4cout << "     daughter 2:" << daughtermomentum[2]/GeV << "[GeV/c]" << G4endl;
    G4cout << "   momentum sum:" << momentumsum/GeV << "[GeV/c]" << G4endl;
  }
#endif

  G4double costheta, sintheta, phi, sinphi, cosphi;
  G4double costhetan, sinthetan, phin, sinphin, cosphin;
  costheta = 2.*G4UniformRand() - 1.0;
  sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta));
  phi  = twopi*G4UniformRand()*rad;
  sinphi = std::sin(phi);
  cosphi = std::cos(phi);
  G4ThreeVector direction0(sintheta*cosphi,sintheta*sinphi,costheta);
  G4DynamicParticle * daughterparticle = new G4DynamicParticle(G4MT_daughters[0], direction0*daughtermomentum[0]);
  products->PushProducts(daughterparticle);

  costhetan = (daughtermomentum[1]*daughtermomentum[1] - daughtermomentum[2]*daughtermomentum[2] - daughtermomentum[0]*daughtermomentum[0])
                                                    /(2.0*daughtermomentum[2]*daughtermomentum[0]);
  sinthetan = std::sqrt((1.0 - costhetan)*(1.0 + costhetan));
  phin  = twopi*G4UniformRand()*rad;
  sinphin = std::sin(phin);
  cosphin = std::cos(phin);
  G4ThreeVector direction2;
  direction2.setX(sinthetan*cosphin*costheta*cosphi - sinthetan*sinphin*sinphi + costhetan*sintheta*cosphi);
  direction2.setY(sinthetan*cosphin*costheta*sinphi + sinthetan*sinphin*cosphi + costhetan*sintheta*sinphi);
  direction2.setZ(-sinthetan*cosphin*sintheta + costhetan*costheta);
  daughterparticle = new G4DynamicParticle(G4MT_daughters[2], direction2*(daughtermomentum[2]/direction2.mag()));
  products->PushProducts(daughterparticle);

  daughterparticle = new G4DynamicParticle(G4MT_daughters[1],
                                           (direction0*daughtermomentum[0] + direction2*(daughtermomentum[2]/direction2.mag()))*(-1.0));
  products->PushProducts(daughterparticle);

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
     G4cout << "SinsqDecayChannel::ThreeBodyDecayIt create decay products in rest frame " << G4endl;
     products->DumpInfo();
  }
#endif
  return products;
}

G4DecayProducts *SinsqDecayChannel::ManyBodyDecayIt()
{
  G4int index, index2;
  G4DecayProducts *products;

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1)
    G4cout << "SinsqDecayChannel::ManyBodyDecayIt()" << G4endl;
#endif

  G4double parentmass = G4MT_parent_mass;
  G4double *daughtermass = new G4double[numberOfDaughters];
  G4double sumofdaughtermass = 0.0;

  for (index=0; index<numberOfDaughters; index++) {
    daughtermass[index] = G4MT_daughters_mass[index];
    sumofdaughtermass += daughtermass[index];
  }
  G4double *daughtermomentum = new G4double[numberOfDaughters];
  G4ThreeVector direction;
  G4DynamicParticle **daughterparticle;
  G4double *sm = new G4double[numberOfDaughters];
  G4double tmas;
  G4double weight = 1.0;
  G4int numberOfTry = 0;

  do {
    G4double temp;
    G4double *rd = new G4double[numberOfDaughters];
    rd[0] = 1.0;

    for (index=1; index<numberOfDaughters-1; index++)
      rd[index] = G4UniformRand();
    rd[numberOfDaughters-1] = 0.0;

    for (index=1; index<numberOfDaughters-1; index++) {
      for (index2=index+1; index2<numberOfDaughters; index2++) {
        if (rd[index] < rd[index2]) {
          temp = rd[index];
          rd[index] = rd[index2];
          rd[index2] = temp;
        }
      }
    }
    tmas = parentmass - sumofdaughtermass;
    temp = sumofdaughtermass;
    for (index=0; index<numberOfDaughters; index++) {
      sm[index] = rd[index]*tmas + temp;
      temp -= daughtermass[index];
      if (GetVerboseLevel() > 1) {
        G4cout << index << "  rundom number:" << rd[index] << "   virtual mass:" << sm[index]/GeV << "[GeV/c/c]" << G4endl;
      }
    }
    delete [] rd;

    weight = 1.0;
    index = numberOfDaughters - 1;
    daughtermomentum[index] = Pmx(sm[index-1], daughtermass[index-1], sm[index]);
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1) {
      G4cout << "     daughter " << index << ":" << *daughters_name[index] << " momentum:" << daughtermomentum[index]/GeV << "[GeV/c]" << G4endl;
    }
#endif
    for (index=numberOfDaughters-2; index>=0; index--) {
      daughtermomentum[index] = Pmx(sm[index],daughtermass[index], sm[index+1]);

      if (daughtermomentum[index] < 0.0) {
#ifdef G4VERBOSE
        if (GetVerboseLevel()>0) {
          G4cout << "SinsqDecayChannel::ManyBodyDecayIt " << "     can not calculate daughter momentum " << G4endl;
          G4cout << "     parent:" << *parent_name << " mass:" << parentmass/GeV << "[GeV/c/c]" << G4endl;
          G4cout << "     daughter " << index << ":" << *daughters_name[index] << " mass:" << daughtermass[index]/GeV << "[GeV/c/c]";
          G4cout << " mass:" << daughtermomentum[index]/GeV << "[GeV/c]" <<G4endl;
        }
#endif
        delete [] sm;
        delete [] daughtermass;
        delete [] daughtermomentum;

        G4Exception("SinsqDecayChannel::ManyBodyDecayIt", "can not create decay products", JustWarning, "sum of daughter mass is larger than parent mass");
        return 0;

      } else {
        weight *= daughtermomentum[index]/sm[index];
#ifdef G4VERBOSE
        if (GetVerboseLevel() > 1) {
          G4cout << "     daughter " << index << ":" << *daughters_name[index] << " momentum:" << daughtermomentum[index]/GeV << "[GeV/c]" << G4endl;
        }
#endif
      }
    }

#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1) {
      G4cout << "    weight: " << weight << G4endl;
    }
#endif

    if (numberOfTry++ > 100) {
#ifdef G4VERBOSE
      if (GetVerboseLevel() > 0) {
        G4cout << "SinsqDecayChannel::ManyBodyDecayIt: " << " can not determine Decay Kinematics " << G4endl << "parent : " << *parent_name << G4endl;
        G4cout << "daughters : ";
        for (index=0; index<numberOfDaughters; index++) {
          G4cout << *daughters_name[index] << " , ";
        }
        G4cout << G4endl;
      }
#endif

      G4Exception("SinsqDecayChannel::ManyBodyDecayIt: ", " Cannot decay ", JustWarning, " Decay Kinematics cannot be calculated ");
      delete [] sm;
      delete [] daughtermass;
      delete [] daughtermomentum;

      return 0;
    }
  } while (weight > G4UniformRand());

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "Start calulation of daughters momentum vector " << G4endl;
  }
#endif

  G4double costheta, sintheta, phi;
  G4double beta;
  daughterparticle = new G4DynamicParticle*[numberOfDaughters];

  index = numberOfDaughters - 2;
  costheta = 2.*G4UniformRand() - 1.0;
  sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta));
  phi = twopi*G4UniformRand()*rad;
  direction.setZ(costheta);
  direction.setY(sintheta*std::sin(phi));
  direction.setX(sintheta*std::cos(phi));
  daughterparticle[index] = new G4DynamicParticle(G4MT_daughters[index], direction*daughtermomentum[index]);
  daughterparticle[index+1] = new G4DynamicParticle(G4MT_daughters[index+1], direction*(-1.0*daughtermomentum[index]));

  for (index=numberOfDaughters-3; index>=0; index--) {
    costheta = 2.*G4UniformRand() - 1.0;
    sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta));
    phi = twopi*G4UniformRand()*rad;
    direction.setZ(costheta);
    direction.setY(sintheta*std::sin(phi));
    direction.setX(sintheta*std::cos(phi));

    beta = daughtermomentum[index];
    beta /= std::sqrt( daughtermomentum[index]*daughtermomentum[index] + sm[index+1]*sm[index+1]);

    for (index2=index+1; index2<numberOfDaughters; index2++) {
      G4LorentzVector p4;
      p4 = daughterparticle[index2]->Get4Momentum();
      p4.boost(direction.x()*beta, direction.y()*beta, direction.z()*beta);
      daughterparticle[index2]->Set4Momentum(p4);
    }
    daughterparticle[index]= new G4DynamicParticle(G4MT_daughters[index], direction*(-1.0*daughtermomentum[index]));
  }
  G4DynamicParticle *parentparticle;
  direction.setX(1.0);
  direction.setY(0.0);
  direction.setZ(0.0);
  parentparticle = new G4DynamicParticle(G4MT_parent, direction, 0.0);
  products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  for (index=0; index<numberOfDaughters; index++) {
    products->PushProducts(daughterparticle[index]);
  }

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "SinsqDecayChannel::ManyBodyDecayIt " << "  create decay products in rest frame " << G4endl;
    products->DumpInfo();
  }
#endif

  delete [] daughterparticle;
  delete [] daughtermomentum;
  delete [] daughtermass;
  delete [] sm;

  return products;
}

G4double SinsqDecayChannel::Pmx(G4double e, G4double p1, G4double p2)
{
  G4double ppp = (e + p1 + p2)*(e + p1 - p2)*(e - p1 + p2)*(e - p1 - p2)/(4.0*e*e);

  if (ppp>0)
    return std::sqrt(ppp);
  else
    return -1.;
}
