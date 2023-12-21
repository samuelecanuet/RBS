// PhysicsList.cc

#include "PhysicsList.hh"
//#include "PhysicsListMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4StepLimiter.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
//#include "G4EmPenelopePhysics.hh"
//#include "G4DecayPhysics.hh"
#include "G4LossTableManager.hh"

#include "G4ProcessManager.hh"
#include "G4PhysicsConstructorRegistry.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronElasticProcess.hh"
#include "G4Proton.hh"
#include "G4EmDNAPhysics.hh"
G4VPhysicsConstructor* GetPhysicsConstructor(const G4String& name) {
	return G4PhysicsConstructorRegistry::Instance()->GetPhysicsConstructor(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() :
	G4VModularPhysicsList() {

	defaultCutValue = 1000 * nm;


	//fMessenger = 0;
	//SetVerboseLevel(1);

	// EM physics
	// RegisterPhysics(GetPhysicsConstructor("G4EmStandard_opt4"));
	// RegisterPhysics(new G4EmStandardPhysics_option4(0));
	// RegisterPhysics( new G4EmDNAPhysics());
	// RegisterPhysics( new G4EmPenelopePhysics());
 	RegisterPhysics( new G4EmLivermorePhysics(1));
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList() {

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name) {
	RegisterPhysics(GetPhysicsConstructor(name));
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
