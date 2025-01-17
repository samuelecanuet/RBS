
#include <iostream>
#include "TFile.h"
#include <fstream>
#include <dirent.h>
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TH1D.h"
#include "TF1.h"
#include <iomanip>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>

using namespace std;
using namespace ROOT::Math;


int THREAD;
bool flag_saving = false;

const vector<string> Face = {"A", "B"};
vector<double> Energy = {1.2, 3.0};
vector<string> Type = {"THIN", "THICK"};
map<string, string> FaceOffSet;


map<string, map<double, map<string, string>>> Exp_FileName_MAP;
map<string, map<double, double>> fMIN_MAP;
map<string, map<double, double>> fMAX_MAP;

map<string, map<double, double>> Sim_WinMIN_MAP;
map<string, map<double, double>> Sim_WinMAX_MAP;
map<string, map<double, double>> Exp_WinMIN_MAP;
map<string, map<double, double>> Exp_WinMAX_MAP;

map<string, double> Thickness_Al1_MAP;
map<string, double> Thickness_Mylar_MAP;
map<string, double> Thickness_Al2_MAP;

map<double, double> Calibration_Coefficient_MAP;
map<double, double> Calibration_Offset_MAP;

map<string, map<double, map<string, double>>> ScaleAl1_MAP;
map<string, map<double, map<string, double>>> ScaleOxygen_MAP;
map<string, map<double, map<string, double>>> ScaleCarbon_MAP;
map<string, map<double, map<string, double>>> ScaleAl2_MAP;

map<string, map<double, map<string, TH1D *>>> Exp_Hist_MAP;

map<string, map<double, map<string, TH1D *>>> Exp_Hist_calib_MAP;
map<string, map<double, map<string, TH1D *>>> Sim_Hist_MAP;
map<string, map<double, map<string, TH1D *>>> Sim_Hist_conv_MAP;

map<string, map<double, map<string, TH1D *>>> Al1_MAP;
map<string, map<double, map<string, TH1D *>>> Oxygen_MAP;
map<string, map<double, map<string, TH1D *>>> Carbon_MAP;
map<string, map<double, map<string, TH1D *>>> Al2_MAP;

map<string, map<double, map<string, double>>> CHI2_MAP;

map<string, double> Optimum_Al;
map<string, double> Optimum_Mylar;

TH1D *Exp_Hist = new TH1D("Exp_Hist", "Experimental RBS", 1024, 0, 1024);
vector<TH1D *> Exp_Hist_AB;
TH1D *Al1;
vector<TH1D *> Al1_AB;
TH1D *Oxygen;
vector<TH1D*> Oxygen_AB;
TH1D *Carbon;
vector<TH1D*> Carbon_AB;
TH1D *Al2;
vector<TH1D*> Al2_AB;

TH1D *Exp_Hist_calib;
vector<TH1D*> Exp_Hist_calib_AB;
TH1D *Sim_Hist;
vector<TH1D*> Sim_Hist_AB;
TH1D *Sim_Hist_conv;
vector<TH1D*> Sim_Hist_conv_AB;

const double fMIN = -1111;
const double fMAX = -1111;

double energy;
TFile *Sim_File;
TFile *fSave;
TFile *f;
string Macro_File;

string MacroModifier(const string& filename, const string& type, double energy, const string& face)
{
    double T1 = Thickness_Al1_MAP[type];
    double T2 = Thickness_Mylar_MAP[type];
    double T3 = Thickness_Al2_MAP[type];

    ostringstream oss;
    oss << fixed << setprecision(1) << energy;
    string new_filename = "../tmp/shoot_" + oss.str() + "_" + to_string((int)T1) + "_" + to_string((int)T2) + "_" + to_string((int)T3) + "_" + face + "_" + FaceOffSet[face] + ".mac";

    ifstream file;
    file.open("../shoot_base.mac");
    ofstream temp_file(new_filename);
    string line;

    
    if (face == "B") 
    {
        swap(T1, T3);
    }

    while (getline(file, line))
    {
        if (line.find("%e") != string::npos)
        {
            line.replace(line.find("%e"), 2, to_string(energy) + " MeV");
        }
        if (line.find("%1") != string::npos)
        {
            line.replace(line.find("%1"), 2, to_string((int)T1) + " nm");
        }
        if (line.find("%2") != string::npos)
        {
            line.replace(line.find("%2"), 2, to_string((int)T2) + " nm");
        }
        if (line.find("%3") != string::npos)
        {
            line.replace(line.find("%3"), 2, to_string((int)T3) + " nm");
        }
        if (line.find("%p") != string::npos)
        {
            line.replace(line.find("%p"), 2, FaceOffSet[face] + " mm");
        }
        temp_file << line << endl;
    }

    file.close();
    temp_file.close();

    return new_filename;
}

string LaunchG4(string macro_filename)
{
    string filename = macro_filename.substr(13, macro_filename.size() - 17) + ".root";
    const string filename_f = "../../../../../../../mnt/hgfs/shared-2/ROOT_files/" + filename;
    system(("cd ../build && example " + macro_filename + " " + filename_f).c_str());

    return filename_f;
}

void CreateRootHist(string Run_Type = "RBS")
{
    TFile *expfile = new TFile("Exp_Hist.root", "RECREATE");
    string data_nb = "15";
    if (Run_Type != "RBS")
    {
        data_nb = "14";
    }
    for (string type : Type)
    {
        for (double energy : Energy)
        {
            for (string face : Face)
            {
                Exp_Hist_MAP[type][energy][face] = new TH1D(("Exp_Hist_" + type + "_" + to_string(energy) + "_" + face).c_str(), ("Experimental " + Run_Type).c_str(), 1024, 0, 1024);
                
                // read txt file
                ifstream file;
                file.open(Exp_FileName_MAP[type][energy][face]);
                string line;
                bool data = false;
                int counter = 0;

                while (getline(file, line))
                {

                    int value;
                    stringstream ss(line);
                    if (line.compare(0, 14, ("[DATA" + data_nb + ",1024 ]").c_str()) == 0)
                    {
                        data = true;
                        continue;
                    }
                    if (data)
                    {
                        ss >> value;
                        Exp_Hist_MAP[type][energy][face]->SetBinContent(counter, value);
                        counter++;
                    }
                    if (line.compare(0, 1, "[") == 0)
                    {
                        data = false;
                    }
                }
                file.close();

                Exp_Hist_MAP[type][energy][face]->Write();
            }
        }
    }
    expfile->Close();
}

void ReadHist()
{
    TFile *f = new TFile("Exp_Hist.root", "READ");
    for (string type : Type)
    {
        for (double energy : Energy)
        {
            for (string face : Face)
            {
                Exp_Hist_MAP[type][energy][face] = (TH1D *)f->Get(("Exp_Hist_" + type + "_" + to_string(energy) + "_" + face).c_str());
            }
        }
    }
}

double Chi2Test(TH1D* h1, TH1D* h2, double fMin, double fMax)
{
    double chi2 = 0;
    for (int i = 0; i < h1->GetNbinsX(); i++)
    {       
        if (h1->GetBinCenter(i) < fMin || h1->GetBinCenter(i) > fMax)
            continue;
        double error = h1->GetBinError(i);
        if (error == h1->GetBinContent(i) + h2->GetBinContent(i))
            continue;
        double diff = h1->GetBinContent(i) - h2->GetBinContent(i);
        chi2 += diff * diff / (h1->GetBinContent(i) + h2->GetBinContent(i));
    }
    return chi2;
}

double FunctionToMinimize(const double *par)
{
    
    Calibration_Offset_MAP[1.2] = par[0];
    Calibration_Coefficient_MAP[1.2] = par[1];
    Calibration_Offset_MAP[3.0] = par[2];
    Calibration_Coefficient_MAP[3.0] = par[3];

    Thickness_Al1_MAP["THIN"] = par[4];
    Thickness_Mylar_MAP["THIN"] = par[5];
    Thickness_Al2_MAP["THIN"] = par[4];

    Thickness_Al1_MAP["THICK"] = par[7];
    Thickness_Mylar_MAP["THICK"] = par[8];
    Thickness_Al2_MAP["THICK"] = par[7];

    double CHI2_SUM = 0;
    double coef = 1;

    for (const auto& type : Type)
    {
        double TH_coef_Al = (int)Thickness_Al1_MAP[type]/Optimum_Al[type];
        double TH_coef_Mylar = (int)Thickness_Mylar_MAP[type]/Optimum_Mylar[type];
        for (const auto& energy : Energy)
        {
            for (const auto& face : Face)
            {
                //////////////////// ############### SIMULATION ############### ////////////////////
                ostringstream osss;
                osss << fixed << setprecision(1) << energy;
                string sss = osss.str();

                ostringstream os_al1;
                os_al1 << fixed << setprecision(1) << Thickness_Al1_MAP[type];

                ostringstream os_mylar;
                os_mylar << fixed << setprecision(1) << Thickness_Mylar_MAP[type];

                Sim_Hist_conv_MAP[type][energy][face] = (TH1D *)fSave->Get(("Sim_Hist_conv_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face] + "_4keV").c_str());
                int integral = Exp_Hist_MAP[type][energy][face]->Integral();
                
                if (!Sim_Hist_conv_MAP[type][energy][face])
                {

                    coef = 1;
                    if (face == "B") // coef from experimental normalization
                    {
                        coef = (double)Exp_Hist_MAP[type][energy]["B"]->Integral() / (double)Exp_Hist_MAP[type][energy]["A"]->Integral();
                    }

                    ostringstream oss;
                    oss << fixed << setprecision(1) << energy;
                    string ss = oss.str();
                    string filename = "../../../../../../../mnt/hgfs/shared-2/ROOT_files/" + oss.str() + "_" + to_string((int)Thickness_Al1_MAP[type]) + "_" + to_string((int)Thickness_Mylar_MAP[type]) + "_" + to_string((int)Thickness_Al2_MAP[type]) + "_" + face + "_" + FaceOffSet[face] + ".root";
                    Sim_File = new TFile(filename.c_str(), "READ");
                    if (Sim_File->IsZombie())
                    {
                        // ofstream file;
                        // file.open("tmp/" + to_string(THREAD) + ".txt");
                        // file << "0";
                        // file.close();
                        // fSave->Close();

                        // return 0;

                        string new_macro_filename = MacroModifier(Macro_File, type, energy, face);
                        LaunchG4(new_macro_filename);
                        Sim_File = new TFile(filename.c_str(), "READ");
                    }

                    Al1 = (TH1D *)Sim_File->Get("RBS_0_13");
                    Oxygen = (TH1D *)Sim_File->Get("RBS_1_8");
                    Carbon = (TH1D *)Sim_File->Get("RBS_1_6");
                    Al2 = (TH1D *)Sim_File->Get("RBS_2_13");

                    Al1->Scale(ScaleAl1_MAP[type][energy][face] / Al1->Integral() * coef * TH_coef_Al);
                    Oxygen->Scale(ScaleOxygen_MAP[type][energy][face] / Oxygen->Integral() * coef * TH_coef_Mylar);
                    Carbon->Scale(ScaleCarbon_MAP[type][energy][face] / Carbon->Integral() * coef * TH_coef_Mylar);
                    Al2->Scale(ScaleAl2_MAP[type][energy][face] / Al2->Integral() * coef * TH_coef_Al);

                    Sim_Hist = (TH1D *)Al1->Clone("Sim_Hist");
                    Sim_Hist->Reset();
                    Sim_Hist->Add(Al1);
                    Sim_Hist->Add(Oxygen);
                    Sim_Hist->Add(Carbon);
                    Sim_Hist->Add(Al2);

                    // convolution sim_hist * gaussian
                    TF1 *Gaussian = new TF1("Gaussian", "gausn", -1000, 1000);
                    Gaussian->SetParameters(1, 0, 4);
                    Sim_Hist_conv = (TH1D *)Sim_Hist->Clone("Sim_Hist_conv");
                    Sim_Hist_conv->Reset();

                    double res;
                    double err;
                    for (int bin = 0; bin < Sim_Hist->GetNbinsX(); bin++)
                    {
                        res = 0;
                        err = 0;
                        for (int bin_c = 0; bin_c < Sim_Hist->GetNbinsX(); bin_c++)
                        {
                            res += Gaussian->Eval(Sim_Hist->GetBinCenter(bin) - Sim_Hist->GetBinCenter(bin_c)) * Sim_Hist->GetBinContent(bin_c);
                            err += Gaussian->Eval(Sim_Hist->GetBinCenter(bin) - Sim_Hist->GetBinCenter(bin_c)) * Sim_Hist->GetBinError(bin_c);
                        }
                        Sim_Hist_conv->SetBinContent(bin, res);
                        Sim_Hist_conv->SetBinError(bin, err);
                    }

                    Al1_MAP[type][energy][face] = (TH1D *)Al1->Clone(("Al1_" + face).c_str());
                    Oxygen_MAP[type][energy][face] = (TH1D *)Oxygen->Clone(("Oxygen_" + face).c_str());
                    Carbon_MAP[type][energy][face] = (TH1D *)Carbon->Clone(("Carbon_" + face).c_str());
                    Al2_MAP[type][energy][face] = (TH1D *)Al2->Clone(("Al2_" + face).c_str());
                    // Sim_Hist_MAP[type][energy][face] = (TH1D *)Sim_Hist->Clone(("Sim_Hist_" + face).c_str());
                    Sim_Hist_conv_MAP[type][energy][face] = (TH1D *)Sim_Hist_conv->Clone(("Sim_Hist_conv_" + type + "_" + oss.str() + "_" + face + "_4keV").c_str());

                    fSave->cd();
                    // Sim_Hist_MAP[type][energy][face]->SetName(("Sim_Hist_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face]).c_str());
                    // Sim_Hist_MAP[type][energy][face]->Write();
                    Sim_Hist_conv_MAP[type][energy][face]->SetName(("Sim_Hist_conv_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face] + "_4keV").c_str());
                    Sim_Hist_conv_MAP[type][energy][face]->Write();
                    // Al1_MAP[type][energy][face]->SetName(("Al1_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face]).c_str());
                    // Al1_MAP[type][energy][face]->Write();
                    // Oxygen_MAP[type][energy][face]->SetName(("Oxygen_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face]).c_str());
                    // Oxygen_MAP[type][energy][face]->Write();
                    // Carbon_MAP[type][energy][face]->SetName(("Carbon_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face]).c_str());
                    // Carbon_MAP[type][energy][face]->Write();
                    // Al2_MAP[type][energy][face]->SetName(("Al2_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face]).c_str());
                    // Al2_MAP[type][energy][face]->Write();

                    // if (!flag_saving)
                    //     Sim_File->Close();
                }
                else
                {
                    // Al1_MAP[type][energy][face] = (TH1D *)fSave->Get(("Al1_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face]).c_str());
                    // Oxygen_MAP[type][energy][face] = (TH1D *)fSave->Get(("Oxygen_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face]).c_str());
                    // Carbon_MAP[type][energy][face] = (TH1D *)fSave->Get(("Carbon_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face]).c_str());
                    // Al2_MAP[type][energy][face] = (TH1D *)fSave->Get(("Al2_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face]).c_str());
                    // Sim_Hist_MAP[type][energy][face] = (TH1D *)fSave->Get(("Sim_Hist_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face]).c_str());
                }
                

                //////////////////// ############### EXPERIMENTAL ############### ////////////////////
                Exp_Hist_calib = (TH1D *)Sim_Hist_conv_MAP[type][energy][face]->Clone("Exp_Hist_calib");
                Exp_Hist_calib->Reset();
                
                for (int bin = 0; bin < integral; bin++)
                {
                    Exp_Hist_calib->Fill(Calibration_Offset_MAP[energy] + Calibration_Coefficient_MAP[energy] * Exp_Hist_MAP[type][energy][face]->GetRandom());
                }
                Exp_Hist_calib_MAP[type][energy][face] = (TH1D *)Exp_Hist_calib->Clone(("Exp_Hist_" + type + "_" + osss.str() + "_" + face).c_str());
                //////////////////// ############### chi2 ############### ////////////////////

                Exp_Hist_calib_MAP[type][energy][face]->GetXaxis()->SetRangeUser(fMIN_MAP[type][energy], fMAX_MAP[type][energy]);
                Sim_Hist_conv_MAP[type][energy][face]->GetXaxis()->SetRangeUser(fMIN_MAP[type][energy], fMAX_MAP[type][energy]);

                // CHI2_MAP[type][energy][face] = Exp_Hist_calib_MAP[type][energy][face]->Chi2Test(Sim_Hist_conv_MAP[type][energy][face], " UU CHI2/NDF");
                CHI2_MAP[type][energy][face] = Chi2Test(Exp_Hist_calib_MAP[type][energy][face], Sim_Hist_conv_MAP[type][energy][face], fMIN_MAP[type][energy], fMAX_MAP[type][energy]) / (fMAX_MAP[type][energy] - fMIN_MAP[type][energy]);

                CHI2_SUM += CHI2_MAP[type][energy][face];
            }
        }
    }

    return CHI2_SUM;
}

int main(int argc, char* argv[])
{
    // converting char to double parameters 
    const double par[10] = {atof(argv[5]), atof(argv[6]), atof(argv[7]), atof(argv[8]), atof(argv[1]), atof(argv[2]), atof(argv[1]), atof(argv[3]), atof(argv[4]), atof(argv[3])};
    int thread = atoi(argv[9]);

    if (par[4] == -1111) {Type.erase(Type.begin());} // thin
    if (par[7] == -1111) {Type.erase(Type.begin() + 1);} // thick

    if (par[1] == -1111) {Energy.erase(Energy.begin());} // 1.2
    if (par[3] == -1111) {Energy.erase(Energy.begin()+1);} // 3.0

    THREAD = thread;
    fSave = new TFile("RBS_Saving.root", "UPDATE");
    f = new TFile("RBS_Results.root", "RECREATE");

    /// EXPERIMENTAL FILE FROM AIFIRA ///
    Exp_FileName_MAP["THICK"][3.0]["A"] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_A_Pos2_3MeV.mpa";
    Exp_FileName_MAP["THICK"][3.0]["B"] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_B_Pos2_3MeV.mpa";
    fMIN_MAP["THICK"][3.0] = 1900;
    fMAX_MAP["THICK"][3.0] = 2700;

    Exp_FileName_MAP["THICK"][1.2]["A"] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_A_Pos2_1,2MeV.mpa";
    Exp_FileName_MAP["THICK"][1.2]["B"] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_B_Pos2_1,2MeV.mpa";
    fMIN_MAP["THICK"][1.2] = 200;
    fMAX_MAP["THICK"][1.2] = 1200;

    Exp_FileName_MAP["THIN"][3.0]["A"] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_A_Pos1_3MeV.mpa";
    Exp_FileName_MAP["THIN"][3.0]["B"] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_B_Pos1_3MeV.mpa";
    fMIN_MAP["THIN"][3.0] = 2100;
    fMAX_MAP["THIN"][3.0] = 2800;

    Exp_FileName_MAP["THIN"][1.2]["A"] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_A_Pos1_1,2MeV.mpa";
    Exp_FileName_MAP["THIN"][1.2]["B"] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_B_Pos1_1,2MeV.mpa";
    fMIN_MAP["THIN"][1.2] = 800;
    fMAX_MAP["THIN"][1.2] = 1100;
    //////////////////////////////////////

    Macro_File = "../shoot.mac";

    // CreateRootHist();
    ReadHist();

    // INITIAL GUESS //

    //////////////////
    Optimum_Al["THICK"] = 110;
    Optimum_Mylar["THICK"] = 6100;

    Optimum_Al["THIN"] = 85;
    Optimum_Mylar["THIN"] = 500;   
    //////////////////


    // OFFSET
    FaceOffSet["A"] = "-2";
    FaceOffSet["B"] = "5";  

    // THICKNESS
    Thickness_Al1_MAP["THIN"] = 100;
    Thickness_Mylar_MAP["THIN"] = 535;
    Thickness_Al2_MAP["THIN"] = 100;

    Thickness_Al1_MAP["THICK"] = 120;
    Thickness_Mylar_MAP["THICK"] = 6090;
    Thickness_Al2_MAP["THICK"] = 120;

    // CALIBRATION
    // Calibration_Offset_MAP[1.2] = 1.89;
    // Calibration_Coefficient_MAP[1.2] = 2.00495;

    // Calibration_Offset_MAP[3.0] = 25;
    // Calibration_Coefficient_MAP[3.0] = 3.7989;

    //moi
    Calibration_Offset_MAP[1.2] = 2.8087;
    Calibration_Coefficient_MAP[1.2] = 2.0032;

    Calibration_Offset_MAP[3.0] = 29.69;
    Calibration_Coefficient_MAP[3.0] = 3.7960;

    // SCALE    
    //// A 
    ScaleAl1_MAP["THIN"][1.2]["A"] = 5000;
    ScaleOxygen_MAP["THIN"][1.2]["A"] = 7400;
    ScaleCarbon_MAP["THIN"][1.2]["A"] = 24500;
    ScaleAl2_MAP["THIN"][1.2]["A"] = 5000;

    ScaleAl1_MAP["THIN"][3.0]["A"] = 2000;
    ScaleOxygen_MAP["THIN"][3.0]["A"] = 6000;
    ScaleCarbon_MAP["THIN"][3.0]["A"] = 15500;
    ScaleAl2_MAP["THIN"][3.0]["A"] = 1000;

    ScaleAl1_MAP["THICK"][1.2]["A"] = 9000;
    ScaleOxygen_MAP["THICK"][1.2]["A"] = 80000;
    ScaleCarbon_MAP["THICK"][1.2]["A"] = 300000;
    ScaleAl2_MAP["THICK"][1.2]["A"] = 15000;

    ScaleAl1_MAP["THICK"][3.0]["A"] = 3700;
    ScaleOxygen_MAP["THICK"][3.0]["A"] = 64000;
    ScaleCarbon_MAP["THICK"][3.0]["A"] = 165000;
    ScaleAl2_MAP["THICK"][3.0]["A"] = 1500;

    //// B
    ScaleAl1_MAP["THIN"][1.2]["B"] = 5000;
    ScaleOxygen_MAP["THIN"][1.2]["B"] = 7200;
    ScaleCarbon_MAP["THIN"][1.2]["B"] = 24000;
    ScaleAl2_MAP["THIN"][1.2]["B"] = 4000;

    ScaleAl1_MAP["THIN"][3.0]["B"] = 2300;
    ScaleOxygen_MAP["THIN"][3.0]["B"] = 6500;
    ScaleCarbon_MAP["THIN"][3.0]["B"] = 15500;
    ScaleAl2_MAP["THIN"][3.0]["B"] = 950;

    ScaleAl1_MAP["THICK"][1.2]["B"] = 7000;
    ScaleOxygen_MAP["THICK"][1.2]["B"] = 80000;
    ScaleCarbon_MAP["THICK"][1.2]["B"] = 300000;
    ScaleAl2_MAP["THICK"][1.2]["B"] = 12000;

    ScaleAl1_MAP["THICK"][3.0]["B"] = 3400;
    ScaleOxygen_MAP["THICK"][3.0]["B"] = 64000;
    ScaleCarbon_MAP["THICK"][3.0]["B"] = 162000;
    ScaleAl2_MAP["THICK"][3.0]["B"] = 2000;
    //////////////////////

    // window
    Sim_WinMIN_MAP["THIN"][1.2] = 1030;
    Sim_WinMAX_MAP["THIN"][1.2] = 1080;
    Exp_WinMIN_MAP["THIN"][1.2] = 512;
    Exp_WinMAX_MAP["THIN"][1.2] = 550;

    Sim_WinMIN_MAP["THIN"][3.0] = 2560;
    Sim_WinMAX_MAP["THIN"][3.0] = 2700;
    Exp_WinMIN_MAP["THIN"][3.0] = 670;
    Exp_WinMAX_MAP["THIN"][3.0] = 710;

    Sim_WinMIN_MAP["THICK"][1.2] = 1020;
    Sim_WinMAX_MAP["THICK"][1.2] = 1080;
    Exp_WinMIN_MAP["THICK"][1.2] = 500;
    Exp_WinMAX_MAP["THICK"][1.2] = 550;

    Sim_WinMIN_MAP["THICK"][3.0] = 2580;
    Sim_WinMAX_MAP["THICK"][3.0] = 2700;
    Exp_WinMIN_MAP["THICK"][3.0] = 660;
    Exp_WinMAX_MAP["THICK"][3.0] = 720;


    // const double par[10] = {
    //     Calibration_Offset_MAP[1.2], Calibration_Coefficient_MAP[1.2],
    //     Calibration_Offset_MAP[3.0], Calibration_Coefficient_MAP[3.0],
    //     Thickness_Al1_MAP["THIN"], Thickness_Mylar_MAP["THIN"], Thickness_Al2_MAP["THIN"],
    //     Thickness_Al1_MAP["THICK"], Thickness_Mylar_MAP["THICK"], Thickness_Al2_MAP["THICK"]
    // };

    flag_saving = true;
    
    double chi2 = FunctionToMinimize(par);
    


    //write chi2 in txt file in thread as name
    ofstream file;
    file.open("tmp/"+to_string(thread) + ".txt");
    for (const auto& type : Type)
    {
        for (const auto& energy : Energy)
        {
            for (const auto& face : Face)
            {
                file << CHI2_MAP[type][energy][face] << endl;
            }
        }
    }
    file.close();

    if (THREAD == 0)
    {
        f->cd();
        for (const auto &type : Type)
        {
            for (const auto &energy : Energy)
            {
                for (const auto &face : Face)
                {
                    Exp_Hist_calib_MAP[type][energy][face]->Write();
                    Sim_Hist_conv_MAP[type][energy][face]->Write();
                }
            }
        }
    }
    f->Close();

    
    fSave->Close();

    return 0;
}