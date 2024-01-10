#include "Layer.hh"
#include "G4VProcess.hh"
#include "G4TrackingManager.hh"
#include "G4Event.hh"
#include <G4Field.hh>
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4LogicalVolume.hh"
#include "G4TouchableHandle.hh"
#include "G4Navigator.hh"

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "globals.hh"


//test fluctuations
#include "G4IonYangFluctuationModel.hh"
#include "G4IonChuFluctuationModel.hh"

#include "G4hIonEffChargeSquare.hh"


//----------------------------------------------------------------------

Layer::Layer(G4int Number, G4Material *Material, G4double Area, G4double Angle, G4double Distance)
    : G4VSensitiveDetector(std::to_string(fNumber)), fNumber(Number), fMaterial(Material), fArea(Area), fAngle(Angle), fDistance(Distance)
{
}

// destructor
Layer::~Layer()
{
}

//----------------------------------------------------------------------

void Layer::Initialize(G4HCofThisEvent *)
{
}

G4bool Layer::ProcessHits(G4Step *step, G4TouchableHistory *)
{
  // G4Track *track = step->GetTrack();
  // if (track->GetParentID() == 0)
  // {
  //   if (step->IsFirstStepInVolume())
  //   {
  //     energy.clear();
  //     position.clear();
  //     direction.clear();
  //   }

  //   if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessType() == fTransportation)
  //   {
  //     InitialParticle = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition();
  //     G4VParticleChange *particleChange = G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager()->GetfParticleChange();
  //     particleChange->SetNumberOfSecondaries(energy.size());
  //     G4double sum = 0;
  //     for (int i = 0.; i < energy.size(); i++)
  //     {
  //       G4ThreeVector Direction = DirectionRBS();
  //       const G4Element *Element = RandomElement(fMaterial);
  //       G4Track *newtrack = new G4Track(new G4DynamicParticle(InitialParticle, Direction, CollisionEnergy(energy[i], Element, Direction.theta())), 0, position[i]);

  //       G4double weight = 1 / (fMaterial->GetDensity() / g * cm * cm * cm) * CalculateTotalRBSYield(energy[i], InitialParticle->GetAtomicMass(), Element->GetA() / (g / mole), InitialParticle->GetAtomicNumber(), Element->GetZ(), Direction.theta());

  //       size_t processCount = InitialParticle->GetProcessManager()->GetProcessList()->size();
  //       G4ProcessVector *vec = InitialParticle->GetProcessManager()->GetProcessList();
  //       G4double crossel = 0;
  //       for (size_t j = 0; j < processCount; ++j)
  //       {
  //         G4VProcess *currentProcess = (*vec)[j];
  //         for (int el = 0; el < fMaterial->GetNumberOfElements(); el++)
  //         {
  //           crossel += fMaterial->GetFractionVector()[el] * emCalculator.ComputeCrossSectionPerAtom(energy[i], "proton", currentProcess->GetProcessName(), fMaterial->GetElement(el), 0.1 * CLHEP::eV) / barn;
  //         }
  //       }

  //       // particleChange->AddSecondary(newtrack);
  //       // particleChange->GetSecondary(i)->SetWeight(crossel/weight);
  //       //G4cout<<crossel<<"        "<<energy[i]<<G4endl;
  //       fRunAction->FillH1(fMaterial->GetName(), position[i].z()/nm);
  //     }
  //   }

  //   if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessType() != fTransportation)
  //   {
      
  //     energy.push_back(step->GetPostStepPoint()->GetKineticEnergy());
  //     position.push_back(step->GetPostStepPoint()->GetPosition());
  //     direction.push_back(step->GetPostStepPoint()->GetMomentumDirection());
  //     step_lenght.push_back(step->GetStepLength());

      
  //   }
  // }

  // else{
  //       step->GetTrack()->SetTrackStatus(fStopAndKill);

  // }
  return (true);
}

G4int Layer::GetNumber()
{
  return fNumber;
}


G4double Layer::CollisionEnergy(G4double E, const G4Element *Element, G4double angle)
{
  // kinematic factor
  //  M1 - incident particle, M2 - target atom, E - incident energy, angle - scattering angle

  G4double M2 = Element->GetAtomicMassAmu();
  G4double M1 = 1;
  G4double k;
  G4double square = std::sqrt(std::pow(M2 / M1, 2.) - std::pow(sin(angle), 2.));
  G4double M1co = cos(angle);
  G4double denominator = 1 + (M2 / M1);
  k = std::pow((M1co + square) / denominator, 2.);
  return E * k;
}

G4ThreeVector Layer::DirectionRBS()
{
  G4ThreeVector dir;
	G4double THETA = THETA = GetRandomTheta(0);
	G4double phi = static_cast<double>(rand()) / RAND_MAX * (-M_PI / 4) - 2.7*pi / 7;

	// G4cout << THETA << "   " << phi*CLHEP::deg
	dir[0] = sin(THETA) * cos(phi);
	dir[1] = sin(THETA) * sin(phi);
	dir[2] = cos(THETA);

	return dir;
}

G4double Layer::GetRandomTheta(G4double alpha)
{
    G4double theta;
    do {
      
        // Générer un nombre aléatoire uniforme dans la plage [minTheta, maxTheta]
        G4double minTheta = (135 - 7 - 2*5) * CLHEP::deg;
        G4double maxTheta = (135  +7+5*2)* CLHEP::deg;
        theta = G4UniformRand() * (maxTheta - minTheta) + minTheta;

        // Calculer le terme sin^4(theta)
        G4double sin4Theta = std::pow(std::sin(theta), -4);

        // Générer un nombre aléatoire y uniformément entre 0 et 1
        G4double y = G4UniformRand()*( std::pow(std::sin(maxTheta), -4) - std::pow(std::sin(minTheta), -4)  ) + std::pow(std::sin(minTheta), -4);

        // Accepter theta si y est inférieur à 1/sin^4(theta)
        if (y <= sin4Theta) {
          
            return theta;
        }
        // Répéter si y > 1/sin^4(theta)
    } while (true);
}


const G4Element *Layer::RandomElement(G4Material *Material)
{
  G4double randomNum = G4UniformRand();
  G4double sum = 0;
  G4double A = 0;
  for (int el = 0; el < Material->GetNumberOfElements(); el++)
  {
    sum += Material->GetFractionVector()[el];
    if (randomNum <= sum)
    {
      return Material->GetElement(el);
    }
  }
  G4Exception("Layer::RandomElement", "ElementNotFound", JustWarning, ("Random element not found in G4Material: " + Material->GetName()).c_str());
  return Material->GetElement(0);
}

////////////////////////////////////////////////////////###################

G4double Layer::CalcAngleCMFrame(G4double angle, G4double M1, G4double M2)
{
  // old version

  G4double CMang = 0., x = 0., y = 0.;
  x = (M1 / M2) * sin(angle);
  y = asin(x);
  CMang = y + angle;

  return CMang;
}

G4double Layer::CalcEnergyCMFrame(G4double energy, G4double M1, G4double M2)
{
  G4double en = 0.;
  en = energy * M2 / (M1 + M2);
  return en;
}

G4double Layer::CalcDiffRuthXsecCM(G4double E, G4double angleCM, G4double Z1, G4double Z2)
{
  G4double sine = sin(angleCM / 2);
  G4double epsilon = 55.26349406 / 1000; // in units of e^2/(MeV*fm)
  G4double FirstTermDenom = 16 * pi * epsilon * E;
  G4double FirstTerm = std::pow((Z1 * Z2) / (FirstTermDenom), 2.);
  G4double SecondTerm = 1 / (std::pow(sine, 4.));
  G4double fullTerm = FirstTerm * SecondTerm;
  return fullTerm * 10;
}

G4double Layer::CalcAndersenScreening(G4double energy_cm, G4double angle_cm, G4double Z1, G4double Z2)
{
  G4double screening_factor;
  G4double V1 = (0.04873 * Z1 * Z2 * std::sqrt(std::pow(Z1, 2. / 3.) + std::pow(Z2, 2. / 3.))) * keV; // in keV
  G4double first_term = std::pow(1 + (0.5 * (V1 / energy_cm)), 2.);
  G4double second_term = V1 / (2 * energy_cm * sin(angle_cm / 2));
  G4double second_term_sq = std::pow(second_term, 2);
  G4double denominator = std::pow(1 + (V1 / energy_cm) + second_term_sq, 2);
  screening_factor = first_term / denominator;

  return screening_factor;
}

G4double Layer::CalcDiffRuthXsecLAB(G4double M1, G4double M2, G4double angle, G4double xsection)
{
  G4double ratio = M1 / M2;
  G4double multiplier = std::pow(1 + (ratio * ratio) + (2 * ratio * cos(angle)), 3. / 2.) / (1 + (ratio * cos(angle)));
  // G4cout << " multiplier " << multiplier << G4endl;
  G4double labXSEC = multiplier * xsection;
  return labXSEC;
}

G4double Layer::CalculateTotalRBSYield(G4double energy, G4double M1, G4double M2, G4double Z1, G4double Z2, G4double angle)
{
  G4double CMangleX = CalcAngleCMFrame(angle, M1, M2);
  G4double CMenergyX = CalcEnergyCMFrame(energy, M1, M2);
  G4double CMxsecX = CalcDiffRuthXsecCM(CMenergyX, CMangleX, Z1, Z2);
  G4double andersen_factor = CalcAndersenScreening(CMenergyX, CMangleX, Z1, Z2);
  G4double CMxsecModX = CMxsecX * andersen_factor;
  G4double CMtoLABxsecX = CalcDiffRuthXsecLAB(M1, M2, CMangleX, CMxsecModX);

  return CMtoLABxsecX;
}


// nuclear energy loss straggling
G4double Layer::CalcNuclEnStraggling(G4double Z1, G4double Z2, G4double M1, G4double M2, G4double atdens,G4double distance)
{
	G4double nucl_strag = 0.;
	nucl_strag = 0.26*std::pow(Z1*Z2,2)*std::pow(M1/(M1+M2),2)*(distance/cm)*(atdens)/(1e+24); // in MeV2

	return nucl_strag;
}

	//energy straggling in MeV^2
G4double Layer::CalcBohrStrag(G4double Z1, G4double Z2, G4double atomDens, G4double dist)
{
	G4double elmcharge_squared = 1.4399764;	// in units of MeV*fm
	G4double e2 = std::pow(elmcharge_squared,2.)*1e-26;	// in units of MeV^2*cm^2
	G4double fourpi = 4*3.14159265358979;
	G4double fBohr = fourpi * e2;
	G4double Bohr_t = Z1*Z1*Z2*fBohr*atomDens*(dist/cm);
	return Bohr_t;
}

// yang+chu
G4double Layer::CalculateTotalBohrStraggling(G4double energy, G4ParticleDefinition* particle, G4Material* mat, G4double distance)
{

	G4double straggling = 0;
	G4double nucl_strag = 0.;
	G4double elec_strag = 0.;
	
	G4int elements    	= mat->GetNumberOfElements();
    G4double *Z2 		= new G4double[elements];
    G4double *M2 		= new G4double[elements];
    G4double *aDensity 	= new G4double[elements];	

	if (energy > 0) {

	    //functions for en loss straggling evaluation
    	G4hIonEffChargeSquare* eff_charge 	= new G4hIonEffChargeSquare(""); // gets effective charge square of ion based on energy
    	G4IonYangFluctuationModel* model_yang	= new G4IonYangFluctuationModel("");
    	G4IonChuFluctuationModel* model_chu 	= new G4IonChuFluctuationModel("");

    	G4double q_squared	= 0.;
    	G4double yang_sq	= 0.;
		G4double chu_sq	= 0.;

		// Particle parameters
		q_squared 		= eff_charge->TheValue(particle,mat,energy);		// effective charge squared
		G4double Z1 		= particle->GetAtomicNumber();
		G4double M1 		= particle->GetAtomicMass();
		G4double q_sqrt 	= sqrt(q_squared);
		G4double q_sqrt_z 	= q_sqrt/Z1;
		
		yang_sq 		= model_yang->TheValue(particle,mat,energy);
		chu_sq			= model_chu->TheValue(particle,mat,energy);

		G4double strag_factor 	= q_sqrt_z*chu_sq+yang_sq;

    	const G4double *atomDensVector = mat->GetVecNbOfAtomsPerVolume();

    	for (int i=0; i<elements; i++) {
    		Z2[i] = mat->GetElement(i)->GetZ();
    		M2[i] = mat->GetElement(i)->GetA()/(g/mole);
    		aDensity[i] = atomDensVector[i]/(1/cm3);
    	}

		for (int j=0; j<elements;j++) {
			elec_strag += CalcBohrStrag(Z1,Z2[j],aDensity[j],distance)*strag_factor;	// electronic straggling straggling, in MeV2
			nucl_strag += CalcNuclEnStraggling(Z1,Z2[j],M1,M2[j],aDensity[j],distance);	// nuclear straggling, in MeV2

			if (elec_strag < 0) 
				elec_strag = 0.;
			if (nucl_strag < 0) 
				nucl_strag = 0.;
		}
	}
	straggling = elec_strag+nucl_strag;

    delete []Z2;
    delete []M2;
    delete []aDensity;

	return straggling;
}