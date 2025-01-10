
#include "G4VDiscreteProcess.hh"
#include "G4VProcess.hh"
#include "RunAction.hh"
#include "G4EmCalculator.hh"
#include "G4Event.hh"

using namespace CLHEP;

class Stepper : public G4VDiscreteProcess {
public:
    Stepper(const G4String& processName, double);
    virtual ~Stepper();

    G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize, G4ForceCondition* condition);
    G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
    // G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&);
    // G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&);
    // G4double AlongStepGetPhysicalInteractionLength(const G4Track&, G4double, G4double, G4double&, G4GPILSelection*);
    // G4double PostStepGetPhysicalInteractionLength(const G4Track&, G4double, G4ForceCondition*);
    // G4double AtRestGetPhysicalInteractionLength(const G4Track&, G4ForceCondition*);

    RunAction* fRunAction;
    double LastHit;

     G4double CollisionEnergy(G4double Energy, const G4Element *Element, G4double Angle);
  G4ThreeVector DirectionRBS();
  G4double CrossSectionRBS(G4double Energy, const G4Element *Element, G4double theta);
  const G4Element *RandomElement(G4Material *Material);

    G4EmCalculator emCalculator;
  G4double GetRandomTheta(G4double a);

  G4double CalcAngleCMFrame(G4double angle, G4double M1, G4double M2);
  G4double CalcEnergyCMFrame(G4double energy, G4double M1, G4double M2);
  G4double CalcDiffRuthXsecLAB(G4double M1, G4double M2, G4double angle, G4double xsection);
  G4double CalcDiffRuthXsecCM(G4double E, G4double angleCM, G4double Z1, G4double Z2);
  G4double CalcDiffRuthXsec(G4double E, G4double M1, G4double M2, G4double angle, G4double Z1, G4double Z2);
  G4double CalculateTotalRBSYield(G4double energy, G4double M1, G4double M2, G4double Z1, G4double Z2, G4double angle);
  G4double CalcAndersenScreening(G4double energy_cm, G4double angle_cm, G4double Z1, G4double Z2);

  G4double GetElementFractionInMaterial(const G4Element* element, const G4Material* material);
  G4double GetElementDensityInMaterial(const G4Element* element, const G4Material* material);



  G4int savedeventID = -1;
  G4double value;
  G4bool Flag = false;
  G4Track *newtrack = nullptr;

};