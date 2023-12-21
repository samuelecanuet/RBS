
/// \file example.cc
/// \brief Main program of the example

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

// #ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
// #else
#include "G4RunManager.hh"
// #endif

#include "G4UImanager.hh"
// #include "QBBC.hh"
#include "PhysicsList.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "Randomize.hh"
#include "G4Threading.hh"
#include "TFile.h"
#include "TTree.h"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv)
{

  G4UIExecutive *ui = 0;
  if (argc == 1)
  {
    ui = new G4UIExecutive(argc, argv);
  }

  std::string filename = "test.root";
  std::string root = "../Result/";
  if (argc == 3)
  {
    filename = root+argv[2];
  }

  G4RunManager *runManager = new G4RunManager;
  // runManager->SetNumberOfThreads(10);

  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new PhysicsList());    
  runManager->SetUserInitialization(new ActionInitialization(filename));

  G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  G4UImanager *UImanager = G4UImanager::GetUIpointer();

  if (!ui)
  {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }
  else
  {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // // Création du fichier final
  //   TFile* finalFile = TFile::Open(filename.c_str(), "UPDATE");

  //   if (!finalFile) {
  //       // Gestion de l'erreur si le fichier final n'est pas ouvert correctement
  //       std::cerr << "Erreur lors de l'ouverture du fichier final." << std::endl;
  //       return 1;
  //   }

  // std::string command = "hadd -f "+filename.substr(0, filename.find_last_of('.'))+"_merged.root"+" "+filename.substr(0, filename.find_last_of('.'))+"_t*.root";
  // std::system(command.c_str());
  // command = "rm "+filename.substr(0, filename.find_last_of('.'))+"_t*.root " + filename;
  // std::system(command.c_str());

  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
