
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "Messenger.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "Tracking.hh"
#include "G4SystemOfUnits.hh"
#include <vector>
#include <cmath>
#include "G4CrossSectionDataStore.hh"
#include "G4NistManager.hh"
#include "G4HadronicProcessStore.hh"

#include "G4VProcess.hh"
#include "G4NistManager.hh"
#include "G4ProcessVector.hh"
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction *eventAction, RunAction *runAction, Tracking *tracking, PrimaryGeneratorAction *generator)
	: G4UserSteppingAction(), fEventAction(eventAction), fRunAction(runAction), fTracking(tracking), fGenerator(generator)
{

	for (int i = 0; i < Hist_sin4->GetNbinsX(); ++i)
	{
		double x = Hist_sin4->GetBinCenter(i + 1);
		double funcValue = 1.0 / std::pow(std::sin(x / 2), 4);
		Hist_sin4->SetBinContent(i + 1, funcValue);
	}

	InitialParticle = G4ParticleTable::GetParticleTable()->FindParticle("proton");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step *step)
{
	
	G4Track *track = step->GetTrack();
	G4TouchableHandle pretouch = step->GetPreStepPoint()->GetTouchableHandle();
	G4TouchableHandle posttouch = step->GetPostStepPoint()->GetTouchableHandle();
	G4double CopyNrB = posttouch->GetCopyNumber();
	G4double CopyNrA = pretouch->GetCopyNumber();

	if (pretouch->GetVolume()->GetLogicalVolume()->GetName() == "World" && track->GetTrackStatus() == fAlive)
	{
		if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() == "RBS" && track->GetParticleDefinition() == InitialParticle)
		{
			Parent = track->GetTrackID();
			G4double Ekin = step->GetPreStepPoint()->GetKineticEnergy() / keV;
			G4double Ekin_vertex = track->GetVertexKineticEnergy() / keV;
			G4double Ekin_res;

			const G4ParticleDefinition *particle = track->GetParticleDefinition();
			G4double weight = track->GetWeight();

			size_t processCount = particle->GetProcessManager()->GetProcessList()->size();
			G4double crossSection = 0;
			G4double crossel = 0;
			G4ProcessVector *vec = particle->GetProcessManager()->GetProcessList();

			std::vector<Data> ParticleVector = fEventAction->GetParticleVector();
			G4double energybeforehit=0;
			for (int i = 0; i < ParticleVector.size(); i++)
			{
				if (step->GetTrack()->GetVertexPosition() == ParticleVector[i].Position)
				{
					energybeforehit = ParticleVector[i].KineticEnergy;
				} 
			}
			
			for (size_t i = 0; i < processCount; ++i)
			{
				G4VProcess *currentProcess = (*vec)[i]; // Accédez directement à l'élément de la liste
				crossSection += emCalculator.GetCrossSectionPerVolume(Ekin_vertex, particle, currentProcess->GetProcessName(), track->GetLogicalVolumeAtVertex()->GetMaterial());

				for (int el = 0; el < track->GetLogicalVolumeAtVertex()->GetMaterial()->GetNumberOfElements(); el++)
				{
					crossel += track->GetLogicalVolumeAtVertex()->GetMaterial()->GetFractionVector()[el] * emCalculator.ComputeCrossSectionPerAtom(energybeforehit, "proton", currentProcess->GetProcessName(), track->GetLogicalVolumeAtVertex()->GetMaterial()->GetElement(el), 0.1*CLHEP::eV)/barn;
				}
			}
			


			// G4cout<<crossSection<<"      "<<weight<<G4endl;

			// if (crossSection && weight)
			// {
			Ekin_res = Ekin + G4RandGauss::shoot(0, 6.21);

			fRunAction->FillNtupleDColumn(0, 0, Ekin_vertex);
			fRunAction->FillNtupleDColumn(0, 1, Ekin);
			fRunAction->FillNtupleDColumn(0, 2, Ekin_res);
			fRunAction->FillNtupleDColumn(0, 3, weight);
			fRunAction->FillNtupleDColumn(0, 4, crossel);
			fRunAction->FillNtupleDColumn(0, 5, track->GetVertexPosition().z() / nm);

			fRunAction->FillH1(0, Ekin_res);
			fRunAction->FillH1(1, Ekin_res, 1 / crossel);
			fRunAction->FillH1(2, Ekin_res,   crossel / weight);
			fRunAction->FillH1(3, Ekin_res, 1/crossel * weight);
			fRunAction->FillH1(4, Ekin_res, weight);
			fRunAction->FillH1(5, track->GetVertexPosition().z() / nm, 1 / crossel);


			fRunAction->FillH1(track->GetLogicalVolumeAtVertex()->GetMaterial()->GetName(), track->GetVertexPosition().z() / nm);

			fRunAction->AddNtupleRow(0);
			//	}
			track->SetTrackStatus(fStopAndKill);
		}
	}
	if (step->GetPostStepPoint()->GetKineticEnergy() == 0)
	{
		if (track->GetParentID() == 0)
		{
			G4VParticleChange *particleChange = G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager()->GetfParticleChange();
			std::vector<Data> ParticleVector = fEventAction->GetParticleVector();
			particleChange->SetNumberOfSecondaries(ParticleVector.size());
			for (int i = 0; i < ParticleVector.size(); i++)
			{
				//  G4cout<<ParticleVector[i].Direction<<"   "<<ParticleVector[i].KineticEnergy<<"   "<<ParticleVector[i].Time<<"   "<<ParticleVector[i].Position<<G4endl;
				//  G4cout << ParticleVector[i].KineticEnergy << "     " << CollisionEnergy(ParticleVector[i].KineticEnergy, ParticleVector[i].Material) << G4endl;
				G4ThreeVector Direction = DirectionRBS();
				G4NistManager *man = G4NistManager::Instance();
				const G4Element *Element = RandomElement(ParticleVector[i].Material);
				G4Track *newtrack = new G4Track(new G4DynamicParticle(InitialParticle, Direction, CollisionEnergy(ParticleVector[i].KineticEnergy, Element, Direction.theta())), ParticleVector[i].Time, ParticleVector[i].Position);
				//G4double info = CrossSectionRBS(ParticleVector[i].KineticEnergy, Element, Direction.theta());


				G4double coef = 1.;
				G4double coef1 = 1;
				for ( int el = 0; el < ParticleVector[i].Material->GetNumberOfElements(); el++)
				{
					const G4Element *ELEMENT = ParticleVector[i].Material->GetElement(el);
					coef += ParticleVector[i].Material->GetFractionVector()[el] * CalculateTotalRBSYield(ParticleVector[i].KineticEnergy, InitialParticle->GetAtomicMass(), ELEMENT->GetA()/(g/mole), InitialParticle->GetAtomicNumber(), ELEMENT->GetZ(), Direction.theta());
					if (ELEMENT = Element)
					{
						coef1 = ParticleVector[i].Material->GetFractionVector()[el];
					}
				}
				
				G4double info = 1/(ParticleVector[i].Material->GetDensity()/g*cm*cm*cm)*CalculateTotalRBSYield(ParticleVector[i].KineticEnergy, InitialParticle->GetAtomicMass(), Element->GetA()/(g/mole), InitialParticle->GetAtomicNumber(), Element->GetZ(), Direction.theta());

				particleChange->AddSecondary(newtrack);
				particleChange->GetSecondary(i)->SetWeight(info);
			}
		}
		track->SetTrackStatus(fStopAndKill);
	}

	if (track->GetParentID() == 0 && step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessType() != fTransportation && CopyNrB > 0)
	{
		fEventAction->AddParticle(step->GetPostStepPoint()->GetPosition(), step->GetPostStepPoint()->GetMomentumDirection(), step->GetPostStepPoint()->GetKineticEnergy(), step->GetPostStepPoint()->GetGlobalTime(), step->GetPostStepPoint()->GetMaterial());
	}

	if (track->GetDefinition()->GetParticleName() == "e-")
	{
		track->SetTrackStatus(fStopAndKill);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// G4double SteppingAction::CollisionEnergy(G4double Energy, const G4Element *Element, G4double theta)
// {
// 	G4double M = Element->GetAtomicMassAmu();

// 	// G4cout<<Energy << "      " << pow( ( 1*cos(theta) + pow(pow(M, 2)-pow(1*sin(theta), 2), 0.5) ) / (1+M) , 2)*Energy << "       " << M << "       "<< G4endl;

// 	return pow((1 * cos(theta) + pow(pow(M, 2) - pow(1 * sin(theta), 2), 0.5)) / (1 + M), 2) * Energy;
// }

G4double SteppingAction::CollisionEnergy(G4double E, const G4Element *Element, G4double angle)
{
	//kinematic factor
	// M1 - incident particle, M2 - target atom, E - incident energy, angle - scattering angle
	G4double M2 = Element->GetAtomicMassAmu();
	G4double M1 = 1;
	G4double k;
	G4double square = std::sqrt(std::pow(M2/M1,2.)-std::pow(sin(angle),2.));
	G4double M1co = cos(angle);
	G4double denominator = 1+(M2/M1);
	k = std::pow((M1co + square)/denominator, 2.);
	return  E*k;
}

G4ThreeVector SteppingAction::DirectionRBS()
{
	G4ThreeVector dir;
	G4double THETA = Hist_sin4->GetRandom();
	G4double phi = static_cast<double>(rand()) / RAND_MAX * (-M_PI / 4) - pi / 3;

	// G4cout << THETA << "   " << phi*CLHEP::deg
	dir[0] = sin(THETA) * cos(phi);
	dir[1] = sin(THETA) * sin(phi);
	dir[2] = cos(THETA);

	return dir;
}

// G4double SteppingAction::CrossSectionRBS(G4double Energy, const G4Element *Element, G4double theta)
// {
// 	G4double Z = Element->GetZ();
// 	G4double A = std::round(Element->GetN());
// 	G4double M = Element->GetAtomicMassAmu();
// 	// G4cout<<pow(Z / Energy, 2) * (pow(1 / Z, 2)) * (1 - (0.049 * pow(Z, 4 / 3)) / Energy)<<G4endl;
// 	return pow(Z / Energy, 2) *( pow(1/sin(theta/2), -4) - 2*pow(1 / M, 2))* (1 - (0.049 * pow(Z, 4 / 3)) / Energy);
// }

G4double SteppingAction::CrossSectionRBS(G4double E, const G4Element *Element, G4double angle)
{
	G4double Z2 = Element->GetZ();
 	G4double M2 = Element->GetAtomicMassAmu();
	G4double M1 = 1;
	G4double Z1 = 1;
	G4double cose = cos(angle);
	G4double sine = sin(angle);
	G4double epsilon = 55.26349406/1000; // in units of e^2/(MeV*fm)
	G4double FirstTermDenom = 8*pi*epsilon*E;
	G4double FirstTerm = std::pow((Z1*Z2)/(FirstTermDenom),2.);
	G4double SecondTerm = 1/(std::pow(sine,4.));
	G4double ThirdTerm = std::pow(M2*cose+std::pow(std::pow(M2,2)-std::pow(M1*sine,2),0.5),2);
	G4double ThirdTermDenom = M2*std::pow(std::pow(M2,2)-std::pow(M1*sine,2),0.5);
	G4double fullTerm = FirstTerm*SecondTerm*ThirdTerm/ThirdTermDenom*10;

	G4double a = 0.049*Z1*std::pow(Z2,4/3);
	G4double b = a/(E/keV);
	G4double c = 1-b;
	G4double result = fullTerm*c;

	return result;//fm^2 conversion to milibarns
}

const G4Element *SteppingAction::RandomElement(G4Material *Material)
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
	G4Exception("SteppingAction::RandomElement", "ElementNotFound", JustWarning, ("Random element not found in G4Material: " + Material->GetName()).c_str());
	return Material->GetElement(0);
}




////////////////////////////////////////////////////////###################

G4double SteppingAction::CalcAngleCMFrame(G4double angle, G4double M1, G4double M2)
{
	// old version

	G4double CMang=0., x= 0., y=0.;
	x = (M1/M2)*sin(angle);
	y = asin(x);
	CMang = y+angle;

	return CMang;
}

G4double SteppingAction::CalcEnergyCMFrame(G4double energy, G4double M1, G4double M2)
{
	G4double en=0.;
	en = energy*M2/(M1+M2);
	return en;
}

G4double SteppingAction::CalcDiffRuthXsecCM(G4double E, G4double angleCM, G4double Z1, G4double Z2)
{
	G4double sine = sin(angleCM/2);
	G4double epsilon = 55.26349406/1000; // in units of e^2/(MeV*fm)
	G4double FirstTermDenom = 16*pi*epsilon*E;
	G4double FirstTerm = std::pow((Z1*Z2)/(FirstTermDenom),2.);
	G4double SecondTerm = 1/(std::pow(sine,4.));
	G4double fullTerm = FirstTerm*SecondTerm;
	return fullTerm*10;
}

G4double SteppingAction::CalcAndersenScreening(G4double energy_cm, G4double angle_cm, G4double Z1, G4double Z2)
{
	G4double screening_factor;
	G4double V1 = (0.04873*Z1*Z2*std::sqrt(std::pow(Z1,2./3.)+std::pow(Z2,2./3.)))*keV;	// in keV
	G4double first_term = std::pow(1+(0.5*(V1/energy_cm)),2.);
	G4double second_term = V1/(2*energy_cm*sin(angle_cm/2));
	G4double second_term_sq = std::pow(second_term,2);
	G4double denominator = std::pow(1+(V1/energy_cm)+second_term_sq,2);
	screening_factor = first_term/denominator;

	return screening_factor;
}


G4double SteppingAction::CalcDiffRuthXsecLAB(G4double M1, G4double M2, G4double angle, G4double xsection)
{
	G4double ratio = M1/M2;
	G4double multiplier = std::pow(1+(ratio*ratio)+(2*ratio*cos(angle)),3./2.)/(1+(ratio*cos(angle)));
	//G4cout << " multiplier " << multiplier << G4endl;
	G4double labXSEC = multiplier * xsection;
	return labXSEC;
}


G4double SteppingAction::CalculateTotalRBSYield(G4double energy, G4double M1, G4double M2, G4double Z1, G4double Z2, G4double angle)
{
	G4double CMangleX 	= CalcAngleCMFrame(angle,M1,M2);
	G4double CMenergyX 	= CalcEnergyCMFrame(energy,M1,M2);
	G4double CMxsecX 	= CalcDiffRuthXsecCM(CMenergyX,CMangleX,Z1,Z2);
	G4double andersen_factor = CalcAndersenScreening(CMenergyX,CMangleX,Z1,Z2);
	G4double CMxsecModX	= CMxsecX*andersen_factor;
	G4double CMtoLABxsecX	= CalcDiffRuthXsecLAB(M1,M2,CMangleX,CMxsecModX);

	return CMtoLABxsecX;

}