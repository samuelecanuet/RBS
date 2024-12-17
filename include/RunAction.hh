#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>
#include <string>
#include "G4AnalysisManager.hh"
#include "TH1D.h"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

using namespace CLHEP;
class G4Run;
class Messenger;
class PrimaryGeneratorAction;
class DetectorConstruction;

class RunAction : public G4UserRunAction
{
  public:
    RunAction(std::string name, PrimaryGeneratorAction* gene, DetectorConstruction* detector);
    virtual ~RunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    std::map<int, int>MAPHIST;
    int counter=0;

    void FillNtupleDColumn(int index_Ntuple, int index, double value){G4AnalysisManager::Instance()->FillNtupleDColumn(index_Ntuple, index, value);}
    void AddNtupleRow(int index){G4AnalysisManager::Instance()->AddNtupleRow(index);}

    void FillH1(int index, double value){G4AnalysisManager::Instance()->FillH1(index, value);}
    void FillH1(int index, double value, double weight){G4AnalysisManager::Instance()->FillH1(index, value, weight);}
    void FillH1(G4String title, double value)
    {
      for (int index = 0; index < G4AnalysisManager::Instance()->GetNofH1s (); index++)
      {
        if (title == G4AnalysisManager::Instance()->GetH1Name(index))
        {
          G4AnalysisManager::Instance()->FillH1(index, value);
        }
      }
    }

    void FillH1(G4String title, double value, double weight)
    {
      for (int index = 0; index < G4AnalysisManager::Instance()->GetNofH1s (); index++)
      {
        if (title == G4AnalysisManager::Instance()->GetH1Name(index))
        {
          G4AnalysisManager::Instance()->FillH1(index, value, weight);
        }
      }
    }
    
    void FillH1Map(int index, double value)
    {
      if (MAPHIST.find(index) != MAPHIST.end())
      {
        G4AnalysisManager::Instance()->FillH1(MAPHIST[index], value);
      }
      else
      {
        G4cout << "Histogram not found" << G4endl;
      }
    }

  G4double GetMaxThickness();

    void SetEnergy(G4double);
    G4double value = 0;

    Messenger* msg;
    G4double MaxThickness;

    PrimaryGeneratorAction* fGene;
    DetectorConstruction* fDetector;


  private:
  std::string filename;

};

