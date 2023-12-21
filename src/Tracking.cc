#include "Tracking.hh"

#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "G4TrackingManager.hh"
#include "G4EventManager.hh"
#include <vector>
#include "TH1D.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

Tracking::Tracking() : G4UserTrackingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

Tracking::~Tracking()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void Tracking::PreUserTrackingAction(const G4Track *)
{
    // Edges.Energy.clear();
    // Edges.Material.clear();
}

void Tracking::PostUserTrackingAction(const G4Track *track)
{
    // if (track->GetParentID() == 0)
    // {
    //     int Z = 0;
    //     int scale = 10;
    //     int bin = (Edges.Energy[0] - Edges.Energy[Edges.Energy.size() - 1]) * scale;
    //     G4cout << bin << G4endl;
    //     double sup = Edges.Energy[0];
    //     double inf = Edges.Energy[Edges.Energy.size() - 1];
    //     TH1D *hist_sigma = new TH1D("hist_sigma", "cross_section", bin, inf, sup);

    //     for (int E = static_cast<int>(std::round((Edges.Energy[Edges.Energy.size() - 1], 1))); E < static_cast<int>(std::round((Edges.Energy[0], 1))); E += 1 / scale)
    //     {
    //         for (int layer = 0; layer < Edges.Energy.size() - 1; layer++)
    //         {
    //             if (Edges.Energy[layer] < E && Edges.Energy[layer + 1] > E)
    //             {
    //                 Z = Edges.Material[layer]->GetZ();
    //             }
    //         }
    //         hist_sigma->Fill(E, CrossSecrionBS(E, Z));
    //     }
    // }
}

// void Tracking::AddEdge(double energy, G4Material *mat)
// {
//     Edges.Energy.push_back(energy);
//     Edges.Material.push_back(mat);
// }

// double Tracking::CrossSecrionBS(double E, double Z)
// {
//     return Z / E;
// }

