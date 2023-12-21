/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4AnalysisManager.hh"
#include "G4RootNtupleFileManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(std::string name)
    : G4UserRunAction()
{
  filename = name;

  msg = new Messenegr();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run *)
{
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
  analysisManager->CreateH1("RBS_detected", "RBS_detected", 1200, 0, 1200);
  analysisManager->CreateH1("RBS_step_corr", "RBS_step_corr", 1200, 0, 1200);
  analysisManager->CreateH1("RBS", "RBS", 1200, 0, 1200);
  analysisManager->CreateH1("RBS2", "RBS2", 1200, 0, 1200);
  analysisManager->CreateH1("RBS3", "RBS3", 1200, 0, 1200);
  analysisManager->CreateH1("z", "z", 1200, 0, 1200);

  G4cout<<msg->GetParticleEnergy()<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run *)
{
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......