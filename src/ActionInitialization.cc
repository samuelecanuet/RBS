/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "Tracking.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization(std::string filename)
 : G4VUserActionInitialization()
{
  name = filename;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
  RunAction* runAction = new RunAction(name);
  SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{

  PrimaryGeneratorAction* generatorAction = new PrimaryGeneratorAction();
  SetUserAction(generatorAction);

  RunAction* runAction = new RunAction(name);
  SetUserAction(runAction);

  Tracking* tracking = new Tracking();
  SetUserAction(tracking);

  EventAction* eventAction = new EventAction(runAction, tracking);
  SetUserAction(eventAction);
  
  SteppingAction* steppingAction = new SteppingAction(eventAction, runAction, tracking, generatorAction);
  SetUserAction(steppingAction);

  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

