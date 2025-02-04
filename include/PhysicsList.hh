//PhysicsList.hh

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"

class G4VPhysicsConstructor;
class G4StepLimiter;
class Messenger;
//class PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VModularPhysicsList
{
public:

  PhysicsList();
  virtual ~PhysicsList();
  //void AddPhysicsList(const G4String& name);
  void ConstructProcess() ;
  void ConstructParticle();

  Messenger *msgg;

  double fAngle;
  void SetAngle(G4double angle) {fAngle = angle;}


private:

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
