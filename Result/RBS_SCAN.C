using namespace std;
using namespace ROOT::Math;

#include <iostream>
#include "TFile.h"
#include <fstream>
#include <dirent.h>
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include <TError.h>

string RED = "\033[1;31m";
string GREEN = "\033[1;32m";
string YELLOW = "\033[1;33m";
string BLUE = "\033[1;34m";
string MAGENTA = "\033[1;35m";
string CYAN = "\033[1;36m";
string WHITE = "\033[1;37m";
string RESET = "\033[0m";

void ProgressCounter(int cEntry, int TotalEntries, string Prefix = "")
{
    cout << BLUE
         << Form(("\r"+Prefix+" : ").c_str())
         << cEntry
         << " / "
         << TotalEntries
         << flush;

  if (cEntry == TotalEntries-1)
  {
    cout << BLUE
         << Form(("\r"+Prefix+" : ").c_str())
         << "Completed "
         << RESET
         << endl;
  }
}

void ProgressBar(int cEntry, int TotalEntries, clock_t start, clock_t Current, string Prefix = "")
{
  if (cEntry % 50 == 0 && cEntry > 1)
  {
    Current = clock();
    const Char_t *Color;
    Double_t Frac = 1.0 * cEntry / TotalEntries;
    Double_t Timeclock = ((double)(Current - start) / CLOCKS_PER_SEC);
    Double_t TimeLeft = Timeclock * (1 / Frac - 1.);
    Color = "\e[1;31m";

    cout << Form(("\r"+Prefix+" Entries : ").c_str())
         << TotalEntries
         << " --- "
         << Form("%4.2f", 100. * cEntry / TotalEntries) << " %"
         << " --- "
         << " Time Left : " << Form("%2d min ", (int)TimeLeft / 60)
         << Form("%02d sec", (int)TimeLeft % 60)
         << flush;
  }

  if (cEntry == TotalEntries-1)
  {
    Current = clock();
    const Char_t *Color;
    Double_t Frac = 1.0 * cEntry / TotalEntries;
    Double_t Timeclock = ((double)(Current - start) / CLOCKS_PER_SEC);
    Double_t TimeLeft = Timeclock * (1 / Frac - 1.);
    Color = "\e[1;31m";
    cout << Form(("\r"+Prefix+" Entries : ").c_str())
         << TotalEntries
         << " --- "
         << Form("%4.2f", 100. * cEntry / TotalEntries) << " %"
         << " --- "
         << " Time Left : " << Form("%2d min ", (int)TimeLeft / 60)
         << Form("%02d sec", (int)TimeLeft % 60)
         << flush;
    cout << endl;
  }
}

bool flag_saving = false;

const vector<string> Face = {"A", "B"};
const vector<double> Energy = {1.2, 3.0};
const vector<string> Type = {"THIN", "THICK"};
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

string MacroModifier(string filename, string type, double energy, string face)
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
        T1 = Thickness_Al2_MAP[type];
        T3 = Thickness_Al1_MAP[type];
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
            }
        }
    }
}

double FunctionToMinimize(const double *par)
{
    
    Calibration_Offset_MAP[1.2] = par[0];
    Calibration_Coefficient_MAP[1.2] = par[1];
    Calibration_Offset_MAP[3.0] = par[2];
    Calibration_Coefficient_MAP[3.0] = par[3];

    // round value to nearest 5 nm
    int step1 = 5;
    int step2 = 5; 
    Thickness_Al1_MAP["THIN"] = round((par[4]) / step1) * step1;
    Thickness_Mylar_MAP["THIN"] = round((par[5]) / step2) * step2;
    Thickness_Al2_MAP["THIN"] = round((par[4]) / step1) * step1;

    Thickness_Al1_MAP["THICK"] = round((par[7]) / step1) * step1;
    Thickness_Mylar_MAP["THICK"] = round((par[8]) / step2) * step2;
    Thickness_Al2_MAP["THICK"] = round((par[7]) / step1) * step1;

    double CHI2_SUM = 0;
    double coef = 1;

    for (const auto& type : Type)
    {
        for (const auto& energy : Energy)
        {
            for (const auto& face : Face)
            {
                int integral = Exp_Hist_MAP[type][energy][face]->Integral();
                //////////////////// ############### SIMULATION ############### ////////////////////
                ostringstream osss;
                osss << fixed << setprecision(1) << energy;
                string sss = osss.str();

                ostringstream os_al1;
                os_al1 << fixed << setprecision(1) << Thickness_Al1_MAP[type];

                ostringstream os_mylar;
                os_mylar << fixed << setprecision(1) << Thickness_Mylar_MAP[type];

                Sim_Hist_conv_MAP[type][energy][face] = (TH1D *)fSave->Get(("Sim_Hist_conv_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face] + "_4keV").c_str());

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

                        return 0;
                        string new_macro_filename = MacroModifier(Macro_File, type, energy, face);
                        LaunchG4(new_macro_filename);
                        Sim_File = new TFile(filename.c_str(), "READ");
                    }

                    Al1 = (TH1D *)Sim_File->Get("RBS_0_13");
                    Oxygen = (TH1D *)Sim_File->Get("RBS_1_8");
                    Carbon = (TH1D *)Sim_File->Get("RBS_1_6");
                    Al2 = (TH1D *)Sim_File->Get("RBS_2_13");

                    Sim_Hist = (TH1D *)Al1->Clone("Sim_Hist");
                    Sim_Hist->Reset();
                    Sim_Hist->Add(Al1, ScaleAl1_MAP[type][energy][face] / Al1->Integral() * coef);
                    Sim_Hist->Add(Oxygen, ScaleOxygen_MAP[type][energy][face] / Oxygen->Integral() * coef);
                    Sim_Hist->Add(Carbon, ScaleCarbon_MAP[type][energy][face] / Carbon->Integral() * coef);
                    Sim_Hist->Add(Al2, ScaleAl2_MAP[type][energy][face] / Al2->Integral() * coef);

                    // convolution sim_hist * gaussian
                    TF1 *Gaussian = new TF1("Gaussian", "gausn", -1000, 1000);
                    Gaussian->SetParameters(1, 0, 4);
                    Sim_Hist_conv = (TH1D *)Sim_Hist->Clone("Sim_Hist_conv");
                    Sim_Hist_conv->Reset();

                    for (int bin = 0; bin < Sim_Hist->GetNbinsX(); bin++)
                    {
                        double res = 0;
                        for (int bin_c = 0; bin_c < Sim_Hist->GetNbinsX(); bin_c++)
                        {
                            res += Gaussian->Eval(Sim_Hist->GetBinCenter(bin) - Sim_Hist->GetBinCenter(bin_c)) * Sim_Hist->GetBinContent(bin_c);
                        }
                        Sim_Hist_conv->SetBinContent(bin, res);
                        Sim_Hist_conv->SetBinError(bin, sqrt(res));
                    }

                    Sim_Hist_conv_MAP[type][energy][face] = (TH1D *)Sim_Hist_conv->Clone(("Sim_Hist_conv_" + face).c_str());

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


                    if (!flag_saving)
                        Sim_File->Close();
                }

                //////////////////// ############### EXPERIMENTAL ############### ////////////////////
                Exp_Hist_calib_MAP[type][energy][face] = (TH1D *)Sim_Hist_conv_MAP[type][energy][face]->Clone("Exp_Hist_calib");
                Exp_Hist_calib_MAP[type][energy][face]->Reset();

                for (int bin = 0; bin < integral; bin++)
                {
                    Exp_Hist_calib_MAP[type][energy][face]->Fill(Calibration_Offset_MAP[energy] + Calibration_Coefficient_MAP[energy] * Exp_Hist_MAP[type][energy][face]->GetRandom());
                }
                //////////////////// ############### chi2 ############### ////////////////////
                Exp_Hist_calib_MAP[type][energy][face]->GetXaxis()->SetRangeUser(fMIN_MAP[type][energy], fMAX_MAP[type][energy]);
                Sim_Hist_conv_MAP[type][energy][face]->GetXaxis()->SetRangeUser(fMIN_MAP[type][energy], fMAX_MAP[type][energy]);

                CHI2_MAP[type][energy][face] = Exp_Hist_calib_MAP[type][energy][face]->Chi2Test(Sim_Hist_conv_MAP[type][energy][face], "CHI2/NDF");
                CHI2_SUM += CHI2_MAP[type][energy][face];
            }
        }
    }

    // PRINTING CHI2
    // cout << "CHI2: " << CHI2_SUM/8 << endl;
    // cout << "Calibration Offset [THIN]: " << Calibration_Offset_MAP[1.2] << " keV" << endl;
    // cout << "Calibration Coefficient [THIN]: " << Calibration_Coefficient_MAP[1.2] << " keV/CH" << endl;
    // cout << "Calibration Offset [THICK]: " << Calibration_Offset_MAP[3.0] << " keV" << endl;
    // cout << "Calibration Coefficient [THICK]: " << Calibration_Coefficient_MAP[3.0] << " keV/CH" << endl;

    // cout << "Al1 Thickness [THIN]: " << Thickness_Al1_MAP["THIN"] << " nm" << endl;
    // cout << "Mylar Thickness [THIN]: " << Thickness_Mylar_MAP["THIN"] << " nm" << endl;
    // cout << "Al2 Thickness [THIN]: " << Thickness_Al2_MAP["THIN"] << " nm" << endl;

    // cout << "Al1 Thickness [THICK]: " << Thickness_Al1_MAP["THICK"] << " nm" << endl;
    // cout << "Mylar Thickness [THICK]: " << Thickness_Mylar_MAP["THICK"] << " nm" << endl;
    // cout << "Al2 Thickness [THICK]: " << Thickness_Al2_MAP["THICK"] << " nm" << endl;
    // for (string type : Type)
    // {
    //     for (double energy : Energy)
    //     {
    //         for (string face : Face)
    //         {
    //             if (CHI2_MAP[type][energy][face] > 10)
    //             {
    //                 cout << "\033[1;31m";
    //             }
    //             cout << "- " << type << "\t" << energy << "\t" << face << "\t" << CHI2_MAP[type][energy][face] << "\033[0m" << endl;
    //         }
    //     }
    // }
    // cout << endl;

    fSave->cd();

    return CHI2_SUM;
}

void RBS_SCAN()
{

    gErrorIgnoreLevel = kWarning;
    gErrorIgnoreLevel = kFatal;

    

    fSave = new TFile("RBS_Saving.root", "UPDATE");

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

    CreateRootHist();

    // INITIAL GUESS //

    // -2 5.5
// CHI2: 1.88214
// Calibration Offset [THIN]: 85.2619 keV
// Calibration Coefficient [THIN]: 1.9697 keV/CH
// Calibration Offset [THICK]: 288.154 keV
// Calibration Coefficient [THICK]: 3.72882 keV/CH
// Al1 Thickness [THIN]: 85 nm
// Mylar Thickness [THIN]: 525 nm
// Al2 Thickness [THIN]: 85 nm
// Al1 Thickness [THICK]: 110 nm
// Mylar Thickness [THICK]: 6100 nm
// Al2 Thickness [THICK]: 110 nm
// - THIN	1.2	A	2.87573
// - THIN	1.2	B	2.10382
// - THIN	3	A	1.34553
// - THIN	3	B	1.49881
// - THICK	1.2	A	2.0028
// - THICK	1.2	B	3.15244
// - THICK	3	A	1.16752
// - THICK	3	B	0.910434


    // OFFSET
    FaceOffSet["A"] = "-2";
    FaceOffSet["B"] = "5.25";  

    // THICKNESS
    Thickness_Al1_MAP["THIN"] = 85;
    Thickness_Mylar_MAP["THIN"] = 525;
    Thickness_Al2_MAP["THIN"] = 85;

    Thickness_Al1_MAP["THICK"] = 110;
    Thickness_Mylar_MAP["THICK"] = 6100;
    Thickness_Al2_MAP["THICK"] = 110;

    //moi
    Calibration_Offset_MAP[1.2] = 85.2619/4;
    Calibration_Coefficient_MAP[1.2] = 1.9697;

    Calibration_Offset_MAP[3.0] = 288.154/4;
    Calibration_Coefficient_MAP[3.0] = 3.72882;

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


    // SCALE PARAMETERS
    int counter = 10;


    vector<double> par_start = {
        -10, 1.82,
        40, 3.57,
        50, 400,
        50, 6000
    };

    vector<double> par_stop = {
        50, 2.12,
        100, 3.87,
        200, 600,
        200, 6200
    };

    const vector<double> par_step = {
        2, 0.01,
        2, 0.01,
        5, 5,
        5, 5
    };

    vector<double> par_optimum = {
        Calibration_Offset_MAP[1.2], Calibration_Coefficient_MAP[1.2],
        Calibration_Offset_MAP[3.0], Calibration_Coefficient_MAP[3.0],
        Thickness_Al1_MAP["THIN"], Thickness_Mylar_MAP["THIN"],
        Thickness_Al1_MAP["THICK"], Thickness_Mylar_MAP["THICK"]
    };



    TH2D *H_Correlation[8][8];
    int TOTAL = 0;

    for (int first = 0; first < 8; first++)
    {
        for (int second = 0; second < 8; second++)
        {
            if (first <= second)
                continue;

            H_Correlation[first][second] = new TH2D(("Correlation_" + to_string(first) + "_" + to_string(second)).c_str(), ("Correlation_" + to_string(first) + "_" + to_string(second)).c_str(), 
            (int)(abs(par_start[first]-par_stop[first])/par_step[first]), par_start[first], par_stop[first], (int)(abs(par_start[second]-par_stop[second])/par_step[second]), par_start[second], par_stop[second]);

            TOTAL += (int)(abs(par_start[first]-par_stop[first])/par_step[first]) * (int)(abs(par_start[second]-par_stop[second])/par_step[second]);
        }
    }

    cout << "TOTAL: " << TOTAL << endl;

    f = new TFile("Correlation.root", "RECREATE");
    // TFile* f_save = new TFile("Correlation_Saved.root", "READ");
    // TTree *tree = (TTree *)f_save->Get("Tree");
    // TTreeReader *Reader = new TTreeReader(tree);
    // TTreeReaderValue<double> al1(Reader, "al1");
    // TTreeReaderValue<double> mylar1(Reader, "mylar1");
    // TTreeReaderValue<double> al2(Reader, "al2");
    // TTreeReaderValue<double> mylar2(Reader, "mylar2");
    // TTreeReaderValue<double> o1(Reader, "o1");
    // TTreeReaderValue<double> c1(Reader, "c1");
    // TTreeReaderValue<double> o3(Reader, "o3");
    // TTreeReaderValue<double> c3(Reader, "c3");
    // TTreeReaderValue<double> chi2(Reader, "chi2");

    // TTreeReaderValue<double> par_saved[8] = {al1, mylar1, al2, mylar2, o1, c1, o3, c3};

    vector<double> par_saving = {0, 0, 0, 0, 0, 0, 0, 0};
    double chi2_saving;
    TFile* f_saving = new TFile("Correlation_Saving.root", "RECREATE");
    TTree* t = new TTree("Tree", "Tree");
    t->Branch("al1", &par_saving[0]);
    t->Branch("mylar1", &par_saving[1]);
    t->Branch("al2", &par_saving[2]);
    t->Branch("mylar2", &par_saving[3]);
    t->Branch("o1", &par_saving[4]);
    t->Branch("c1", &par_saving[5]);
    t->Branch("o3", &par_saving[6]);
    t->Branch("c3", &par_saving[7]);
    t->Branch("chi2", &chi2_saving);

    

    TCanvas *c = new TCanvas("c", "c", 800, 800);
    c->Divide(8, 8);

    clock_t start = clock(), Current;
    int count = 0;
    int counter_progress = 0;
    for (int first = 0; first < 8; first++)
    {
        for (int second = 0; second < 8; second++)
        {
            count++;
            if (first <= second)
                continue;
            
            for (double par1 = par_start[first]; par1 < par_stop[first]; par1 += par_step[first])
            {
                for (double par2 = par_start[second]; par2 < par_stop[second]; par2 += par_step[second])
                {
                    vector<double> par_it = par_optimum;
                    par_it[first] = par1;
                    par_it[second] = par2;

                    const double par[10] = {
                        par_it[0], par_it[1],
                        par_it[2], par_it[3],
                        par_it[4], par_it[5], 0,
                        par_it[6], par_it[7], 0,

                    };

                    flag_saving = false;
                    double Chi2 = FunctionToMinimize(par);

                    H_Correlation[first][second]->Fill(par1, par2, Chi2);

                    chi2_saving = Chi2;
                    par_saving = par_it;
                    t->Fill();
                    t->AutoSave("FlushBaskets");
                    

                    counter_progress++;
                    // ProgressCounter(counter_progress, TOTAL, "Progress: ");
                    ProgressBar(counter_progress, TOTAL, start, Current, "Progress: ");
                }
                f->cd();
                H_Correlation[first][second]->Write("", TObject::kOverwrite);
            }

            c->cd(count);
            H_Correlation[first][second]->Draw("colz");
        }
    }

    c->SaveAs("Correlation.pdf");
    
    f->cd();
    c->Write();
    f->Close();

    fSave->Close();
}




// -2 5.5
// CHI2: 1.88214
// Calibration Offset [THIN]: 85.2619 keV
// Calibration Coefficient [THIN]: 1.9697 keV/CH
// Calibration Offset [THICK]: 288.154 keV
// Calibration Coefficient [THICK]: 3.72882 keV/CH
// Al1 Thickness [THIN]: 85 nm
// Mylar Thickness [THIN]: 525 nm
// Al2 Thickness [THIN]: 85 nm
// Al1 Thickness [THICK]: 110 nm
// Mylar Thickness [THICK]: 6100 nm
// Al2 Thickness [THICK]: 110 nm
// - THIN	1.2	A	2.87573
// - THIN	1.2	B	2.10382
// - THIN	3	A	1.34553
// - THIN	3	B	1.49881
// - THICK	1.2	A	2.0028
// - THICK	1.2	B	3.15244
// - THICK	3	A	1.16752
// - THICK	3	B	0.910434




// -2 5.25
// CHI2: 2.01541
// Calibration Offset [THIN]: 21.3188 keV
// Calibration Coefficient [THIN]: 1.96968 keV/CH
// Calibration Offset [THICK]: 72.0188 keV
// Calibration Coefficient [THICK]: 3.7288 keV/CH
// Al1 Thickness [THIN]: 85 nm
// Mylar Thickness [THIN]: 550 nm
// Al2 Thickness [THIN]: 85 nm
// Al1 Thickness [THICK]: 110 nm
// Mylar Thickness [THICK]: 6100 nm
// Al2 Thickness [THICK]: 110 nm
// - THIN	1.2	A	3.25823
// - THIN	1.2	B	2.54433
// - THIN	3	A	1.38651
// - THIN	3	B	1.57088
// - THICK	1.2	A	2.05241
// - THICK	1.2	B	3.16409
// - THICK	3	A	1.17723
// - THICK	3	B	0.96959

// Errors: 
// 1.04201
// 0.00142783
// 2.69953
// 0.00269693
// 8.68452
// 5.82943
// 0
// 3.90085
// 2.68783
// 0


