
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

SteppingAction::SteppingAction(EventAction *eventAction, RunAction *runAction, Tracking *tracking, PrimaryGeneratorAction* generator)
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
			G4ProcessVector* vec = particle->GetProcessManager()->GetProcessList();
			for (size_t i = 0; i < processCount; ++i) {
				G4VProcess* currentProcess = (*vec)[i]; // Accédez directement à l'élément de la liste
				crossSection += emCalculator.GetCrossSectionPerVolume(Ekin_vertex, particle, currentProcess->GetProcessName(), track->GetLogicalVolumeAtVertex()->GetMaterial());
			}
			//G4cout<<crossSection<<"      "<<weight<<G4endl;

			// if (crossSection && weight)
			// {
				Ekin_res = Ekin + G4RandGauss::shoot(0, 6.21) ;

				fRunAction->FillNtupleDColumn(0, 0, Ekin_vertex);
				fRunAction->FillNtupleDColumn(0, 1, Ekin);
				fRunAction->FillNtupleDColumn(0, 2, Ekin_res);
				fRunAction->FillNtupleDColumn(0, 3, weight);
				fRunAction->FillNtupleDColumn(0, 4, crossSection);
				fRunAction->FillNtupleDColumn(0, 5, track->GetVertexPosition().z() / nm);

				fRunAction->FillH1(0, Ekin_res);
				fRunAction->FillH1(1, Ekin_res, 1/crossSection);
				fRunAction->FillH1(2, Ekin_res, 1/crossSection * weight);
				fRunAction->FillH1(3, Ekin_res, crossSection / weight);
				fRunAction->FillH1(4, Ekin_res, crossSection * weight);
				fRunAction->FillH1(5, track->GetVertexPosition().z() / nm, 1/crossSection);


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
				G4double info = CrossSectionRBS(ParticleVector[i].KineticEnergy, Element);
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

G4double SteppingAction::CollisionEnergy(G4double Energy, const G4Element *Element, G4double theta)
{
	G4double M = Element->GetAtomicMassAmu();

	// G4cout<<Energy << "      " << pow( ( 1*cos(theta) + pow(pow(M, 2)-pow(1*sin(theta), 2), 0.5) ) / (1+M) , 2)*Energy << "       " << M << "       "<< G4endl;

	return pow((1 * cos(theta) + pow(pow(M, 2) - pow(1 * sin(theta), 2), 0.5)) / (1 + M), 2) * Energy;
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

G4double SteppingAction::CrossSectionRBS(G4double Energy, const G4Element *Element)
{
	G4double Z = Element->GetZ();
	G4double A = std::round(Element->GetN());
	// G4cout<<pow(Z / Energy, 2) * (pow(1 / Z, 2)) * (1 - (0.049 * pow(Z, 4 / 3)) / Energy)<<G4endl;
	return pow(Z / Energy, 2) * (1 - (0.049 * pow(Z, 4 / 3)) / Energy);
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

