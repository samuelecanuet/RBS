
/// \file ActionInitialization.hh
/// \brief Definition of the ActionInitialization class

#ifndef ActionInitialization_h
#define ActionInitialization_h 1
#include <string>
#include "G4VUserActionInitialization.hh"
#include "Messenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
using namespace std;
/// Action initialization class.
class ActionInitialization : public G4VUserActionInitialization
{
  public:
    ActionInitialization(std::string filename, DetectorConstruction* detector);
    virtual ~ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
    std::string name;

    PrimaryGeneratorAction* GetGene();

    Messenger* msg;
    DetectorConstruction* fDetector;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
