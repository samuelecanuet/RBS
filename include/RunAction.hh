/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>
#include <string>
#include "G4AnalysisManager.hh"

class G4Run;
class Messenger;

class RunAction : public G4UserRunAction
{
  public:
    RunAction(std::string name);
    virtual ~RunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    std::string GenerateUniqueFilenameByDateTime();

    void FillNtupleDColumn(int index_Ntuple, int index, double value){G4AnalysisManager::Instance()->FillNtupleDColumn(index_Ntuple, index, value);}
    void AddNtupleRow(int index){G4AnalysisManager::Instance()->AddNtupleRow(index);}
    void FillH1(int index, double value){G4AnalysisManager::Instance()->FillH1(index, value);}
    void FillH1(int index, double value, double weight){G4AnalysisManager::Instance()->FillH1(index, value, weight);}


    Messenger* msg;


  private:
  std::string filename;

};

#endif

