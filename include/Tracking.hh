#include "G4UserTrackingAction.hh"
#include "G4Material.hh"
#include <vector>



class Tracking : public G4UserTrackingAction
{
public:
  Tracking();

  virtual ~Tracking();

  virtual void PreUserTrackingAction(const G4Track *);
  virtual void PostUserTrackingAction(const G4Track *);

  // Edges_struct Edges;
  // void AddEdge(double e, G4Material* mat);
  // double CrossSecrionBS(double e, double Z);
};
