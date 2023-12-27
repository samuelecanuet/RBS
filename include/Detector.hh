
#include "G4VSensitiveDetector.hh"
#include "G4EventManager.hh"

class Detector : public G4VSensitiveDetector
{

  //------------------------------------------------------------
public:
  Detector();
  ~Detector();

  void Initialize(G4HCofThisEvent *);
  G4bool ProcessHits(G4Step *, G4TouchableHistory *);
};
