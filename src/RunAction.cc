/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4AnalysisManager.hh"
#include "G4RootNtupleFileManager.hh"
#include "Messenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(std::string name, PrimaryGeneratorAction* gene, DetectorConstruction* detector)
    : G4UserRunAction(), fGene(gene), fDetector(detector)
{
  filename = name; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run *)
{
  
  G4double MaxEnergy = fGene->GetE()/keV;
  std::vector<G4Material*> vec = fDetector->GetMaterials();

  

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(0);

  // Open an output file
  //
  analysisManager->OpenFile(filename);
  analysisManager->CreateNtuple("Catcher", "Hits");
  analysisManager->CreateNtupleDColumn("E_vertex");
  analysisManager->CreateNtupleDColumn("RBS");
  analysisManager->CreateNtupleDColumn("RBS_res");
  analysisManager->CreateNtupleDColumn("crosssection_rbs");
  analysisManager->CreateNtupleDColumn("crosssection_step");
  analysisManager->CreateNtupleDColumn("z_vertex");
  
  analysisManager->FinishNtuple();
  analysisManager->CreateH1("RBS_detected", "RBS_detected", MaxEnergy, 0, MaxEnergy);
  analysisManager->CreateH1("RBS_step_corr", "RBS_step_corr", MaxEnergy, 0, MaxEnergy);
  analysisManager->CreateH1("RBS", "RBS", MaxEnergy, 0, MaxEnergy);
  analysisManager->CreateH1("RBS2", "RBS2", MaxEnergy, 0, MaxEnergy);
  analysisManager->CreateH1("RBS3", "RBS3", MaxEnergy, 0, MaxEnergy);
  analysisManager->CreateH1("z", "z", MaxEnergy, 0, MaxEnergy);

  for (int i=0; i<vec.size(); i++)
  {
    analysisManager->CreateH1(vec[i]->GetName(), vec[i]->GetName(), MaxEnergy, 0, MaxEnergy);
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run *)
{
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
