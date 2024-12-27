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
#include "G4EmStandardPhysicsGS.hh"
//#include "G4EmPenelopePhysics.hh"
//#include "G4DecayPhysics.hh"
#include "G4LossTableManager.hh"

#include "G4ProcessManager.hh"
#include "G4PhysicsConstructorRegistry.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronElasticProcess.hh"
#include "G4Proton.hh"
#include "G4EmSaturation.hh"

#include "Stepper.hh"
G4VPhysicsConstructor* GetPhysicsConstructor(const G4String& name) {
	return G4PhysicsConstructorRegistry::Instance()->GetPhysicsConstructor(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() :
	G4VModularPhysicsList() {

	defaultCutValue = 1. * mm;

	//SetVerboseLevel(1);

	// EM physics
	// RegisterPhysics(GetPhysicsConstructor("G4EmStandard_opt4"));
	// RegisterPhysics(new G4EmStandardPhysics_option4(0));
	// RegisterPhysics( new G4EmDNAPhysics());
	// RegisterPhysics( new G4EmPenelopePhysics());
 	//RegisterPhysics( new G4EmLivermorePhysics(0));

	
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList() {

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void PhysicsList::AddPhysicsList(const G4String& name) {
// 	RegisterPhysics(GetPhysicsConstructor(name));
    
// }

void PhysicsList::ConstructProcess() {
    //Appeler la méthode de la classe de base
    //G4VModularPhysicsList::ConstructProcess();

	AddTransportation();

  // G4VPhysicsConstructor * emPhysicsList = new G4EmPenelopePhysics();
  // G4VPhysicsConstructor * emPhysicsList = new G4EmLivermorePhysics();
  // G4VPhysicsConstructor * emPhysicsList = new G4EmStandardPhysics(1);
  // G4VPhysicsConstructor * emPhysicsList = new G4EmStandardPhysics_option4(0);
  G4VPhysicsConstructor *emPhysicsList = new G4EmStandardPhysicsGS(0);
  emPhysicsList->ConstructProcess();

  G4EmParameters *emParams = G4EmParameters::Instance();
  emParams->SetNumberOfBinsPerDecade(200);

    // Ajouter votre processus personnalisé
    G4ProcessManager* pmanager = G4Proton::Proton()->GetProcessManager();
    if (pmanager) {
        pmanager->AddProcess(new Stepper("CREATOR_PROCESS", fAngle), -1, -1, 1); // Ajouter votre processus personnalisé ici
        // pmanager->SetProcessOrderingToLast(new Stepper("CREATOR_PROCESS"), idxAtRest);
        // pmanager->SetProcessOrderingToLast(new Stepper("CREATOR_PROCESS"), idxAlongStep);
        // pmanager->SetProcessOrderingToLast(new Stepper("CREATOR_PROCESS"), idxPostStep);
    }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{

  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();

  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  // mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();

  // barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();

  // ions
  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4Alpha::AlphaDefinition();
  G4GenericIon::GenericIonDefinition();
}