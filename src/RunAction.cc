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
#include "Detector.hh"
#include "Layer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(std::string name, PrimaryGeneratorAction *gene, DetectorConstruction *detector)
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
  G4double MaxEnergy = fGene->GetE() / keV;
  std::vector<G4Material *> Materials = fDetector->GetMaterials();
  std::vector<Detector *> Detectors = fDetector->GetDetectors();
  std::vector<Layer *> Layers = fDetector->GetLayers();
  MaxThickness = fDetector->GetTotalThickness()/nm;

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(0);

  // Open an output file
  //
  analysisManager->OpenFile(filename);

  // LAYER
  analysisManager->CreateNtuple("Layers", "Hits");
  analysisManager->CreateNtupleDColumn("E_vertex");
  analysisManager->CreateNtupleDColumn("crosssection_rbs");
  analysisManager->CreateNtupleDColumn("crosssection_step");
  analysisManager->CreateNtupleDColumn("z_vertex");
  analysisManager->FinishNtuple(0);


  // DETECTORS
  analysisManager->CreateNtuple("Detectors", "Detectors Hits");
  for (int i = 0; i < Detectors.size(); i++)
  {
    G4String Name = Detectors[i]->GetName();
    analysisManager->CreateNtupleDColumn(Name+"_Ekin");
    analysisManager->CreateNtupleDColumn(Name+"_Ekin_res");
    analysisManager->CreateNtupleDColumn(Name+"_weight");

    analysisManager->CreateH1(Name+"_Ekin_res", Name+"_Ekin_res", MaxEnergy, 0, MaxEnergy);
    analysisManager->CreateH1(Name+"_Ekin_weight", Name+"_Ekin_weight", MaxEnergy, 0, MaxEnergy);
    analysisManager->CreateH1(Name+"_Ekin_1weight", Name+"_Ekin_1weight", MaxEnergy, 0, MaxEnergy);
  }
  analysisManager->FinishNtuple(1);

  for (int i = 0; i < Materials.size(); i++)
  {
    analysisManager->CreateH1(Materials[i]->GetName(), Materials[i]->GetName(), MaxThickness, 0, MaxThickness);
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

G4double RunAction::GetMaxThickness()
{
  G4cout << MaxThickness;
  return MaxThickness;
}
