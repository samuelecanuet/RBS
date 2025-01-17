#include "Stepper.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>
#include "globals.hh"
#include "G4Material.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"
#include "Messenger.hh"
#include "Randomize.hh"
#include "G4EventManager.hh"

Stepper::Stepper(const G4String &processName, double angle) : G4VDiscreteProcess(processName)
{
}

Stepper::~Stepper() {}

G4double Stepper::GetMeanFreePath(const G4Track &aTrack, G4double previousStepSize, G4ForceCondition *condition)
{

    G4double newStepSize = 1 * CLHEP::nm;

    if (aTrack.GetTouchable()->GetVolume()->GetCopyNo() > 0 && aTrack.GetParentID() == 0)
    {
        *condition = Forced;
        return newStepSize;
    }

    return DBL_MAX;
}

G4VParticleChange *Stepper::PostStepDoIt(const G4Track &aTrack, const G4Step &aStep)
{
    G4VParticleChange *particleChange = new G4VParticleChange();

    // if (savedeventID < G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID())
    // {
    //     savedeventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
    //     LastHit = -999;
    // }

    if (!Flag)
    {
        // cout << (int)LastHit << "   " << (int)(aStep.GetPostStepPoint()->GetPosition().z() / nm) << endl;
        if (aStep.GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "CREATOR_PROCESS" && aTrack.GetCurrentStepNumber() < 1e5)
        {
            // LastHit = aStep.GetPostStepPoint()->GetPosition().z()/nm;
            Flag = false;
            particleChange->SetNumberOfSecondaries(1);
            G4Material *fMaterial = aTrack.GetMaterial();
            const G4ParticleDefinition *InitialParticle = aTrack.GetParticleDefinition();
            G4double Energy = aTrack.GetKineticEnergy();

            G4ThreeVector Direction = DirectionRBS();
            const G4Element *Element = RandomElement(fMaterial);
            newtrack = new G4Track(new G4DynamicParticle(InitialParticle, Direction, CollisionEnergy(Energy, Element, Direction.theta())), 0, aTrack.GetPosition());
            // G4cout << "    " << GetElementFractionInMaterial(Element, fMaterial) <<G4endl;

            G4double weight = GetElementDensityInMaterial(Element, fMaterial) / fMaterial->GetDensity() * CalculateTotalRBSYield(Energy, InitialParticle->GetAtomicMass(), Element->GetA() / (g / mole), InitialParticle->GetAtomicNumber(), Element->GetZ(), Direction.theta());

            particleChange->AddSecondary(newtrack);
            // particleChange->GetSecondary(0)->SetWeight(weight);
            particleChange->GetSecondary(0)->SetWeight(aStep.GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo() * 1000 + Element->GetZ());
            // fRunAction->FillH1(fMaterial->GetName(), aTrack.GetPosition().z() / nm);
        }
    }
    return particleChange;
}

G4double Stepper::GetElementFractionInMaterial(const G4Element *element, const G4Material *material)
{
    double fraction = 0.0;

    if (element && material)
    {
        size_t numberOfElements = material->GetNumberOfElements();

        for (size_t i = 0; i < numberOfElements; ++i)
        {
            const G4Element *matElement = material->GetElement(i);
            if (matElement == element)
            {
                fraction = material->GetFractionVector()[i];
                break;
            }
        }
    }

    return fraction;
}

G4double Stepper::GetElementDensityInMaterial(const G4Element *element, const G4Material *material)
{
    double fraction = 0.0;

    if (element && material)
    {
        size_t numberOfElements = material->GetNumberOfElements();
        const G4double *vec1 = material->GetVecNbOfAtomsPerVolume();

        for (size_t i = 0; i < numberOfElements; ++i)
        {
            const G4Element *matElement = material->GetElement(i);
            if (matElement == element)
            {
                fraction = vec1[i] / (1 / cm3);
                break;
            }
        }
    }

    return fraction;
}

G4double Stepper::CollisionEnergy(G4double E, const G4Element *Element, G4double angle)
{
    // kinematic factor
    //  M1 - incident particle, M2 - target atom, E - incident energy, angle - scattering angle

    G4double M2 = Element->GetAtomicMassAmu();
    G4double M1 = 1.00727647;
    G4double k;
    G4double square = std::sqrt(std::pow(M2 / M1, 2.) - std::pow(sin(angle), 2.));
    G4double M1co = cos(angle);
    G4double denominator = 1 + (M2 / M1);
    k = std::pow((M1co + square) / denominator, 2.);
    return E * k;
}

G4ThreeVector Stepper::DirectionRBS()
{
    G4ThreeVector dir;
    G4double THETA = GetRandomTheta(0);

    G4double phi = -90 * CLHEP::deg;
    G4double phiMin = phi - 15 * CLHEP::deg;
    G4double phiMax = phi + 15 * CLHEP::deg;
    phi = phiMin + (phiMax - phiMin) * G4UniformRand();

    // G4cout << THETA << "   " << phi*CLHEP::deg
    dir[0] = sin(THETA) * cos(phi);
    dir[1] = sin(THETA) * sin(phi);
    dir[2] = cos(THETA);

    return dir;
}

G4double Stepper::GetRandomTheta(G4double alpha)
{
    G4double theta = -135;
    G4double minTheta = (theta - 20) * CLHEP::deg;
    G4double maxTheta = (theta + 20) * CLHEP::deg;
    return acos(G4UniformRand() * (cos(maxTheta) - cos(minTheta)) + cos(minTheta));
}

const G4Element *Stepper::RandomElement(G4Material *Material)
{
    G4double randomNum = G4UniformRand();
    G4double sum = 0;
    G4double A = 0;
    G4double num = Material->GetNumberOfElements();
    for (int el = 0; el < Material->GetNumberOfElements(); el++)
    {
        const G4Element *element = Material->GetElement(el);
        sum += element->GetN() * Material->GetFractionVector()[el];
        // G4cout << element->GetN() << "    " << Material->GetAtomicNumDensityVector()[el] << G4endl;
    }
    G4double incr = 0;
    for (int el = 0; el < Material->GetNumberOfElements(); el++)
    {
        const G4Element *element = Material->GetElement(el);

        incr += element->GetN() * Material->GetFractionVector()[el];
        if (randomNum <= incr / sum)
        {
            return Material->GetElement(el);
        }
    }
    // for (int el = 0; el < Material->GetNumberOfElements(); el++)
    // {
    //     sum += 1/num;
    //     if (randomNum <= sum)
    //     {
    //         return Material->GetElement(el);
    //     }
    // }
    G4Exception("Stepper::RandomElement", "ElementNotFound", JustWarning, ("Random element not found in G4Material: " + Material->GetName()).c_str());
    return Material->GetElement(0);
}

////////////////////////////////////////////////////////###################

G4double Stepper::CalcAngleCMFrame(G4double angle, G4double M1, G4double M2)
{
    // old version

    G4double CMang = 0., x = 0., y = 0.;
    x = (M1 / M2) * sin(angle);
    y = asin(x);
    CMang = y + angle;

    return CMang;
}

G4double Stepper::CalcEnergyCMFrame(G4double energy, G4double M1, G4double M2)
{
    G4double en = 0.;
    en = energy * M2 / (M1 + M2);
    return en;
}

G4double Stepper::CalcDiffRuthXsecCM(G4double E, G4double angleCM, G4double Z1, G4double Z2)
{
    G4double sine = sin(angleCM / 2);
    G4double epsilon = 55.26349406 / 1000; // in units of e^2/(MeV*fm)
    G4double FirstTermDenom = 16 * pi * epsilon * E;
    G4double FirstTerm = std::pow((Z1 * Z2) / (FirstTermDenom), 2.);
    G4double SecondTerm = 1 / (std::pow(sine, 4.));
    G4double fullTerm = FirstTerm * SecondTerm;
    return fullTerm * 10;
}

G4double Stepper::CalcDiffRuthXsec(G4double E, G4double M1, G4double M2, G4double angle, G4double Z1, G4double Z2)
{
    G4double cose = cos(angle);
    G4double sine = sin(angle);
    G4double epsilon = 55.26349406 / 1000; // in units of e^2/(MeV*fm)
    G4double FirstTermDenom = 8 * pi * epsilon * E;
    G4double FirstTerm = std::pow((Z1 * Z2) / (FirstTermDenom), 2.);
    G4double SecondTerm = 1 / (std::pow(sin(angle / 2), 4.));
    G4double ThirdTerm = std::pow(M2 * cose + std::pow(std::pow(M2, 2) - std::pow(M1 * sine, 2), 0.5), 2);
    G4double ThirdTermDenom = M2 * std::pow(std::pow(M2, 2) - std::pow(M1 * sine, 2), 0.5);
    G4double fullTerm = FirstTerm * SecondTerm * ThirdTerm / ThirdTermDenom;
    return fullTerm * 10; // fm^2 conversion to milibarns
}

G4double Stepper::CalcAndersenScreening(G4double energy_cm, G4double angle_cm, G4double Z1, G4double Z2)
{
    G4double screening_factor;
    G4double V1 = (0.04873 * Z1 * Z2 * std::sqrt(std::pow(Z1, 2. / 3.) + std::pow(Z2, 2. / 3.))) * keV; // in keV
    G4double first_term = std::pow(1 + (0.5 * (V1 / energy_cm)), 2.);
    G4double second_term = V1 / (2 * energy_cm * sin(angle_cm / 2));
    G4double second_term_sq = std::pow(second_term, 2);
    G4double denominator = std::pow(1 + (V1 / energy_cm) + second_term_sq, 2);
    screening_factor = first_term / denominator;

    return screening_factor;
}

G4double Stepper::CalcDiffRuthXsecLAB(G4double M1, G4double M2, G4double angle, G4double xsection)
{
    G4double ratio = M1 / M2;
    G4double multiplier = std::pow(1 + (ratio * ratio) + (2 * ratio * cos(angle)), 3. / 2.) / (1 + (ratio * cos(angle)));
    // G4cout << " multiplier " << multiplier << G4endl;
    G4double labXSEC = multiplier * xsection;
    return labXSEC;
}

G4double Stepper::CalculateTotalRBSYield(G4double CMenergyX, G4double M1, G4double M2, G4double Z1, G4double Z2, G4double CMangleX)
{
    // G4double CMangleX = CalcAngleCMFrame(angle, M1, M2);
    // G4double CMenergyX = CalcEnergyCMFrame(energy, M1, M2);
    // G4double CMxsecX = CalcDiffRuthXsecCM(CMenergyX, CMangleX, Z1, Z2);
    G4double CMxsecX = CalcDiffRuthXsec(CMenergyX, M1, M2, CMangleX, Z1, Z2);
    G4double andersen_factor = CalcAndersenScreening(CMenergyX, CMangleX, Z1, Z2);
    G4double CMxsecModX = CMxsecX * andersen_factor;
    G4double CMtoLABxsecX = CalcDiffRuthXsecLAB(M1, M2, CMangleX, CMxsecModX);

    return CMxsecX;
}
