
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
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TFitResult.h" 

using namespace std;
using namespace ROOT::Math;


int THREAD;
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


map<string, map<double, map<string, vector<double>>>> Fitting_Params_MAP;
map<string, map<double, map<string, vector<double>>>> Fitting_ParamsError_MAP;
map<double, TGraphErrors*> G_Calibration_MAP;


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

double FunctionFit1(double *x, double *par)
{
    double A1 = par[0];
    double mu11 = par[1];
    double sigma1 = par[2];
    double mu12 = par[3];
    double A2 = par[4];
    double mu21 = par[5];
    double sigma2 = par[6];
    double mu22 = par[7];
    double A3 = par[8];
    double mu3 = par[9];
    double sigma3 = par[10];
    double A4 = par[11];
    double mu4 = par[12];
    double sigma4 = par[13];

    return A1 * (0.5 + erf((x[0] - mu11) / sigma1) / 2) * (erfc((x[0] - mu12) / sigma1) / 2) +
           A2 * (0.5 + erf((x[0] - mu21) / sigma2) / 2) * (erfc((x[0] - mu22) / sigma2) / 2) +
           A3 * exp(-pow((x[0] - mu3) / sigma3, 2)) +
           A4 * exp(-pow((x[0] - mu4) / sigma4, 2));
}

double FunctionFit2(double *x, double *par)
{
    double a = par[0];
    double b = par[1];
    double mu1 = par[2];
    double sigma1 = par[3];
    double mu2 = par[4];
    double sigma2 = par[5];
    double A0 = par[6];
    double mu3 = par[7];
    double sigma3 = par[8];
    double mu4 = par[9];
    double A5 = par[10];
    double mu5 = par[11];
    double sigma5 = par[12];
    double A6 = par[13];
    double mu6 = par[14];
    double sigma6 = par[15];


    return  (a * x[0] + b) * (0.5 + erf((x[0] - mu1) / sigma1) / 2) * (erfc((x[0] - mu2) / sigma2) / 2) +
            A0 * (0.5 + erf((x[0] - mu3) / sigma3) / 2) * (erfc((x[0] - mu4) / sigma3) / 2) +
            A5 * exp(-pow((x[0] - mu5) / sigma5, 2)) +
            A6 * exp(-pow((x[0] - mu6) / sigma6, 2));
}

void FitExp()
{
    TFile *ff = new TFile("RBS_Fit.root", "RECREATE");
    for (string type : Type)
    {
        for (double energy : Energy)
        {
            if (energy == 3.0)
                continue;
            for (string face : Face)
            {
                if (type == "THIN")
                {
                    Exp_Hist = Exp_Hist_MAP[type][energy][face];
                    // Exp_Hist->GetXaxis()->SetRangeUser(fMIN_MAP[type][energy], fMAX_MAP[type][energy]);
                    TF1 *f1 = new TF1("f1", FunctionFit1, 0, 1024, 14);
                    f1->SetNpx(100000);
                    f1->SetParameter(0, 1259);   // A1
                    f1->SetParameter(1, 425);    // mu11
                    f1->SetParameter(2, 6.28);   // sigma1
                    f1->SetParameter(3, 445);    // mu12
                    f1->SetParameter(4, 398.37); // A2
                    f1->SetParameter(5, 458.61); // mu21
                    f1->SetParameter(6, 6.47);   // sigma2
                    f1->SetParameter(7, 479.3);  // mu22
                    f1->SetParameter(8, 590.2);  // A3
                    f1->SetParameter(9, 502);    // mu3
                    f1->SetParameter(10, 3.75);  // sigma3
                    f1->SetParameter(11, 647.1); // A4
                    f1->SetParameter(12, 525);   // mu4
                    f1->SetParameter(13, 3.18);  // sigma4
                    Exp_Hist->Fit("f1", "E");

                    ff->cd();
                    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
                    Exp_Hist->Draw();
                    f1->Draw("same");
                    c1->Write();

                    // cout << "params: " << endl;
                    // cout << "mu11: " << f1->GetParError(1) << endl;
                    // cout << "mu12: " << f1->GetParError(3) << endl;
                    // cout << "mu21: " << f1->GetParError(5) << endl;
                    // cout << "mu22: " << f1->GetParError(7) << endl;
                    // cout << "mu3: " << f1->GetParError(9) << endl;
                    // cout << "mu4: " << f1->GetParError(12) << endl;

                    Fitting_Params_MAP[type][energy][face].push_back(f1->GetParameter(1));
                    Fitting_Params_MAP[type][energy][face].push_back(f1->GetParameter(3));
                    Fitting_Params_MAP[type][energy][face].push_back(f1->GetParameter(5));
                    Fitting_Params_MAP[type][energy][face].push_back(f1->GetParameter(7));
                    Fitting_Params_MAP[type][energy][face].push_back(f1->GetParameter(9));
                    Fitting_Params_MAP[type][energy][face].push_back(f1->GetParameter(12));

                    delete f1;
                }
                else
                {
                    Exp_Hist = Exp_Hist_MAP[type][energy][face];
                    // Exp_Hist->GetXaxis()->SetRangeUser(fMIN_MAP[type][energy], fMAX_MAP[type][energy]);
                    TF1 *f2 = new TF1("f2", FunctionFit2, 0, 1024, 16);
                    f2->SetNpx(100000);
                    f2->SetParameter(0, 0.1);   // a
                    f2->SetParLimits(0, 0, 10);
                    f2->SetParameter(1, 0);   // b
                    
                    f2->SetParameter(2, 200);   // mu1
                    f2->SetParLimits(2, 100, 300);

                    f2->SetParameter(3, 50);  // sigma1
                    f2->SetParLimits(3, 10, 100);

                    f2->SetParameter(4, 445);   // mu2
                    f2->SetParLimits(4, 400, 500);

                    f2->SetParameter(5, 6.47);  // sigma2
                    f2->SetParLimits(5, 1, 50);

                    f2->SetParameter(6, 398.37); // A0
                    f2->SetParLimits(6, 0, 1000);

                    f2->SetParameter(7, 425);   // mu3
                    f2->SetParLimits(7, 400, 450);

                    f2->SetParameter(8, 3.75);  // sigma3
                    f2->SetParLimits(8, 1, 50);

                    f2->SetParameter(9, 525);   // mu4
                    f2->SetParLimits(9, 450, 500);

                    f2->SetParameter(10, 647.1); // A5
                    f2->SetParLimits(10, 0, 1000);

                    f2->SetParameter(11, 250); // mu5
                    f2->SetParLimits(11, 200, 300);

                    f2->SetParameter(12, 3.18);  // sigma5
                    f2->SetParLimits(12, 1, 50);

                    f2->SetParameter(13, 590.2); // A6
                    f2->SetParLimits(13, 0, 1500);

                    f2->SetParameter(14, 510); // mu6
                    f2->SetParLimits(14, 500, 550);

                    f2->SetParameter(15, 3.18);  // sigma6
                    f2->SetParLimits(15, 1, 50);
                    Exp_Hist->Fit("f2", "RE", "", 150, 600);

                    ff->cd();
                    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
                    Exp_Hist->Draw();
                    f2->Draw("same");
                    c1->Write();

                    // cout << "params: " << endl;
                    // cout << "mu1 " << f2->GetParameter(2) << endl;
                    // cout << "mu2 " << f2->GetParameter(4) << endl;
                    // cout << "mu3 " << f2->GetParameter(9) << endl;
                    // cout << "mu5 " << f2->GetParameter(11) << endl;
                    // cout << "mu6 " << f2->GetParameter(14) << endl;


                    // cout << "mu1: " << f2->GetParError(2) << endl;
                    // cout << "mu2: " << f2->GetParError(4) << endl;
                    // cout << "mu3: " << f2->GetParError(9) << endl;
                    // cout << "mu4: " << f2->GetParError(11) << endl;
                    // cout << "mu5: " << f2->GetParError(14) << endl;

                    

                    Fitting_Params_MAP[type][energy][face].push_back(f2->GetParameter(2));
                    Fitting_Params_MAP[type][energy][face].push_back(f2->GetParameter(4));
                    Fitting_Params_MAP[type][energy][face].push_back(f2->GetParameter(9));
                    Fitting_Params_MAP[type][energy][face].push_back(f2->GetParameter(11));
                    Fitting_Params_MAP[type][energy][face].push_back(f2->GetParameter(14));

                    delete f2;
                }
            }
        }
    }
    // ff->Close();    
}

double Chi2Test(TH1D* h1, TH1D* h2, double fMin, double fMax)
{
    double chi2 = 0;
    for (int i = 0; i < h1->GetNbinsX(); i++)
    {       
        if (h1->GetBinCenter(i) < fMin || h1->GetBinCenter(i) > fMax)
            continue;
        double error = h1->GetBinError(i);
        if (error == 0)
            continue;
        double diff = h1->GetBinContent(i) - h2->GetBinContent(i);
        chi2 += diff * diff / error / error;
    }
    return chi2;
}

double FunctionToMinimize(const double *par)
{
    
    Calibration_Offset_MAP[1.2] = par[0];
    // Calibration_Offset_MAP[1.2] = 0;
    Calibration_Coefficient_MAP[1.2] = par[1];
    Calibration_Offset_MAP[3.0] = par[2];
    // Calibration_Offset_MAP[3.0] = 0;
    Calibration_Coefficient_MAP[3.0] = par[3];

    // round value to nearest 5 nm
    const int step1 = 5;
    const int step2 = 5; 
    Thickness_Al1_MAP["THIN"] = round((par[4]) / step1) * step1;
    Thickness_Mylar_MAP["THIN"] = round((par[5]) / step2) * step2;
    // Thickness_Al2_MAP["THIN"] = par[6];
    Thickness_Al2_MAP["THIN"] = round((par[4]) / step1) * step1;

    Thickness_Al1_MAP["THICK"] = round((par[7]) / step1) * step1;
    Thickness_Mylar_MAP["THICK"] = round((par[8]) / step2) * step2;
    // Thickness_Al2_MAP["THICK"] = par[9];
    Thickness_Al2_MAP["THICK"] = round((par[7]) / step1) * step1;

    // ScaleAl1_MAP["THIN"][1.2] = par[10];
    // ScaleOxygen_MAP["THIN"][1.2] = par[11];
    // ScaleCarbon_MAP["THIN"][1.2] = par[12];
    // ScaleAl2_MAP["THIN"][1.2] = par[13];

    // ScaleAl1_MAP["THIN"][3.0] = par[14];
    // ScaleOxygen_MAP["THIN"][3.0] = par[15];
    // ScaleCarbon_MAP["THIN"][3.0] = par[16];
    // ScaleAl2_MAP["THIN"][3.0] = par[17];

    // ScaleAl1_MAP["THICK"][1.2] = par[18];
    // ScaleOxygen_MAP["THICK"][1.2] = par[19];
    // ScaleCarbon_MAP["THICK"][1.2] = par[20];
    // ScaleAl2_MAP["THICK"][1.2] = par[21];

    // ScaleAl1_MAP["THICK"][3.0] = par[22];
    // ScaleOxygen_MAP["THICK"][3.0] = par[23];
    // ScaleCarbon_MAP["THICK"][3.0] = par[24];
    // ScaleAl2_MAP["THICK"][3.0] = par[25];


    double CHI2_SUM = 0;
    double coef = 1;

    G_Calibration_MAP[1.2] = new TGraphErrors();

    for (const auto& type : Type)
    {
        if (type == "THIN")
            continue;
        for (const auto& energy : Energy)
        {
            if (energy == 3.0)
                continue;
            for (const auto& face : Face)
            {
                cout << "type: " << type << " energy: " << energy << " face: " << face << endl;
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

                        ofstream file;
                        file.open("tmp/" + to_string(THREAD) + ".txt");
                        file << "0";
                        file.close();

                        fSave->Close();

                        return 0;

                        string new_macro_filename = MacroModifier(Macro_File, type, energy, face);
                        LaunchG4(new_macro_filename);
                        Sim_File = new TFile(filename.c_str(), "READ");
                    }

                    Al1 = (TH1D *)Sim_File->Get("RBS_0_13");
                    Oxygen = (TH1D *)Sim_File->Get("RBS_1_8");
                    Carbon = (TH1D *)Sim_File->Get("RBS_1_6");
                    Al2 = (TH1D *)Sim_File->Get("RBS_2_13");

                    Al1->Scale(ScaleAl1_MAP[type][energy][face] / Al1->Integral() * coef);
                    Oxygen->Scale(ScaleOxygen_MAP[type][energy][face] / Oxygen->Integral() * coef);
                    Carbon->Scale(ScaleCarbon_MAP[type][energy][face] / Carbon->Integral() * coef);
                    Al2->Scale(ScaleAl2_MAP[type][energy][face] / Al2->Integral() * coef);

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
                    for (int bin = 0; bin < Sim_Hist->GetNbinsX(); bin++)
                    {
                        res = 0;
                        for (int bin_c = 0; bin_c < Sim_Hist->GetNbinsX(); bin_c++)
                        {
                            res += Gaussian->Eval(Sim_Hist->GetBinCenter(bin) - Sim_Hist->GetBinCenter(bin_c)) * Sim_Hist->GetBinContent(bin_c);
                        }
                        Sim_Hist_conv->SetBinContent(bin, res);
                        Sim_Hist_conv->SetBinError(bin, sqrt(res));
                    }

                    Al1_MAP[type][energy][face] = (TH1D *)Al1->Clone(("Al1_" + face).c_str());
                    Oxygen_MAP[type][energy][face] = (TH1D *)Oxygen->Clone(("Oxygen_" + face).c_str());
                    Carbon_MAP[type][energy][face] = (TH1D *)Carbon->Clone(("Carbon_" + face).c_str());
                    Al2_MAP[type][energy][face] = (TH1D *)Al2->Clone(("Al2_" + face).c_str());
                    // Sim_Hist_MAP[type][energy][face] = (TH1D *)Sim_Hist->Clone(("Sim_Hist_" + face).c_str());
                    Sim_Hist_conv_MAP[type][energy][face] = (TH1D *)Sim_Hist_conv->Clone(("Sim_Hist_conv_" + face).c_str());

                    fSave->cd();
                    // Sim_Hist_MAP[type][energy][face]->SetName(("Sim_Hist_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face]).c_str());
                    // Sim_Hist_MAP[type][energy][face]->Write();
                    // Sim_Hist_conv_MAP[type][energy][face]->SetName(("Sim_Hist_conv_" + osss.str() + "_" + os_al1.str() + "_" + os_mylar.str() + "_" + face + "_" + FaceOffSet[face] + "_4keV").c_str());
                    // Sim_Hist_conv_MAP[type][energy][face]->Write();
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

                //////////////////// ############### FITTING ############### ////////////////////
                ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);

                TF1 *f1 = new TF1("f1", FunctionFit1, fMIN_MAP[type][energy], fMAX_MAP[type][energy], 14);
                Sim_Hist_conv_MAP[type][energy][face]->GetXaxis()->SetRangeUser(fMIN_MAP[type][energy], fMAX_MAP[type][energy]);

                if (type == "THIN")
                {
                    // ShowPeaks for guess values, searching 4 peaks.
                    TSpectrum *spectrum = new TSpectrum(4, 1); // Max number of peaks to search for

                    // Find the peaks in the histogram
                    Int_t nPeaks = spectrum->Search(Sim_Hist_conv_MAP[type][energy][face], 10, "", 0.1); // Parameters: hist, sigma, options, threshold

                    // sort the peaks in vector
                    vector<double> peaks;
                    for (int i = 0; i < nPeaks; i++)
                    {
                        // cout << "Peak " << i << " at " << spectrum->GetPositionX()[i] << endl;
                        peaks.push_back(spectrum->GetPositionX()[i]);
                    }
                    sort(peaks.begin(), peaks.end());

                    f1->SetParameter(0, 625); // A1
                    f1->SetParLimits(0, 400, 1000);
                    f1->SetParameter(1, peaks[0] - 20); // mu11
                    f1->SetParLimits(1, peaks[0] - 50, peaks[0] - 10);
                    f1->SetParameter(2, 10.64); // sigma1
                    f1->SetParLimits(2, 2, 20);
                    f1->SetParameter(3, peaks[0] + 20); // mu12
                    f1->SetParLimits(3, peaks[0] + 10, peaks[0] + 50);
                    f1->SetParameter(4, 200); // A2
                    f1->SetParLimits(4, 0, 1000);
                    f1->SetParameter(5, peaks[1] - 15); // mu21
                    f1->SetParLimits(5, peaks[1] - 25, peaks[1] - 5);
                    f1->SetParameter(6, 10.64); // sigma2
                    f1->SetParLimits(6, 2, 20);
                    f1->SetParameter(7, peaks[1] + 15); // mu22
                    f1->SetParLimits(7, peaks[1] + 5, peaks[1] + 25);
                    f1->SetParameter(8, 200); // A3
                    f1->SetParLimits(8, 0, 1000);
                    f1->SetParameter(9, peaks[2]); // mu3
                    f1->SetParLimits(9, peaks[2] - 10, peaks[2] + 10);
                    f1->SetParameter(10, 10.64); // sigma3
                    f1->SetParLimits(10, 2, 20);
                    f1->SetParameter(11, 200); // A4
                    f1->SetParLimits(11, 0, 1000);
                    f1->SetParameter(12, peaks[3]); // mu4
                    f1->SetParLimits(12, peaks[3] - 10, peaks[3] + 10);
                    f1->SetParameter(13, 10.64); // sigma4
                    f1->SetParLimits(13, 2, 20);
                    TFitResultPtr res_ptr = Sim_Hist_conv_MAP[type][energy][face]->Fit("f1", "RES", "", fMIN_MAP[type][energy], fMAX_MAP[type][energy]);

                    cout << res_ptr->IsValid() << endl;

                    // if (res_ptr->IsValid())
                    // {

                        G_Calibration_MAP[energy]->AddPoint(Fitting_Params_MAP[type][energy][face][0], f1->GetParameter(1));
                        // G_Calibration_MAP[energy]->SetPointError(G_Calibration_MAP[energy]->GetN() - 1, sqrt(pow(Fitting_ParamsError_MAP[type][energy][face][0], 2) + pow(Fitting_ParamsError_MAP[type][energy][face][1], 2)), sqrt(pow(f1->GetParError(1), 2) + pow(f1->GetParError(3), 2)));
                        G_Calibration_MAP[energy]->AddPoint(Fitting_Params_MAP[type][energy][face][1], f1->GetParameter(3));
                        // G_Calibration_MAP[energy]->SetPointError(G_Calibration_MAP[energy]->GetN() - 1, Fitting_ParamsError_MAP[type][energy][face][1], f1->GetParError(3));
                        G_Calibration_MAP[energy]->AddPoint(Fitting_Params_MAP[type][energy][face][2], f1->GetParameter(5));
                        // G_Calibration_MAP[energy]->SetPointError(G_Calibration_MAP[energy]->GetN() - 1, sqrt(pow(Fitting_ParamsError_MAP[type][energy][face][2], 2) + pow(Fitting_ParamsError_MAP[type][energy][face][3], 2)), sqrt(pow(f1->GetParError(5), 2) + pow(f1->GetParError(7), 2)));
                        G_Calibration_MAP[energy]->AddPoint(Fitting_Params_MAP[type][energy][face][3], f1->GetParameter(7));
                        // G_Calibration_MAP[energy]->SetPointError(G_Calibration_MAP[energy]->GetN() - 1, Fitting_ParamsError_MAP[type][energy][face][3], f1->GetParError(7));
                        G_Calibration_MAP[energy]->AddPoint(Fitting_Params_MAP[type][energy][face][4], f1->GetParameter(9));
                        // G_Calibration_MAP[energy]->SetPointError(G_Calibration_MAP[energy]->GetN() - 1, Fitting_ParamsError_MAP[type][energy][face][4], f1->GetParError(9));
                        G_Calibration_MAP[energy]->AddPoint(Fitting_Params_MAP[type][energy][face][5], f1->GetParameter(12));
                        // G_Calibration_MAP[energy]->SetPointError(G_Calibration_MAP[energy]->GetN() - 1, Fitting_ParamsError_MAP[type][energy][face][5], f1->GetParError(12));
                    // }
                }
                else
                {
                    TF1 *f2 = new TF1("f2", FunctionFit2, fMIN_MAP[type][energy], fMAX_MAP[type][energy], 16);

                    f2->SetParameter(0, 0.74);     //  a
                    f2->SetParLimits(0, 0., 1);
                    f2->SetParameter(1, 50.583);      //  b
                    f2->SetParLimits(1, 0, 300);
                    f2->SetParameter(2, 400);      //  mu1
                    f2->SetParLimits(2, 300, 500);
                    f2->SetParameter(3, 40.9368);      //  sigma1
                    f2->SetParLimits(3, 10, 80);
                    f2->SetParameter(4, 876.777);      //  mu2
                    f2->SetParLimits(4, 825, 900);
                    f2->SetParameter(5, 15.7857);      //  sigma2
                    f2->SetParLimits(5, 5, 50);
                    f2->SetParameter(6, 150.785);      //  A0
                    f2->SetParLimits(6, 10, 300);
                    f2->SetParameter(7, 900.347);      //  mu3  
                    f2->SetParLimits(7, 850, 940);
                    f2->SetParameter(8, 10.1946);      //  sigma3  
                    f2->SetParLimits(8, 2, 20);
                    f2->SetParameter(9, 960.22);      //   mu4
                    f2->SetParLimits(9, 900, 1000);
                    f2->SetParameter(10, 150.6191);      //  A5
                    f2->SetParLimits(10, 50, 250);
                    f2->SetParameter(11, 560.);      //  mu5
                    f2->SetParLimits(11, 500, 650);
                    f2->SetParameter(12, 33.7196);      //  simga5
                    f2->SetParLimits(12, 2, 60);
                    f2->SetParameter(13, 375.383);     //  A6
                    f2->SetParLimits(13, 200, 1000);
                    f2->SetParameter(14, 1052.79);     //  mu6
                    f2->SetParLimits(14, 1000, 1100);
                    f2->SetParameter(15, 6.89182);     //  sigma6
                    f2->SetParLimits(15, 2, 10);

                    TFitResultPtr res_ptr =  Sim_Hist_conv_MAP[type][energy][face]->Fit("f2", "RES", "", fMIN_MAP[type][energy], fMAX_MAP[type][energy]);
                    // cout << "params: " << endl;
                    // cout << "mu1 " << f2->GetParameter(2) << endl;
                    // cout << "mu2 " << f2->GetParameter(4) << endl;
                    // cout << "mu3 " << f2->GetParameter(9) << endl;
                    // cout << "mu5 " << f2->GetParameter(11) << endl;
                    // cout << "mu6 " << f2->GetParameter(14) << endl;

                    cout << res_ptr->IsValid() << endl;

                    // don't fill the graph if the fit failed
                    if (res_ptr->IsValid())
                    {
                        // G_Calibration_MAP[energy]->AddPoint(Fitting_Params_MAP[type][energy][face][0], f2->GetParameter(2));
                        // G_Calibration_MAP[energy]->SetPointError(G_Calibration_MAP[energy]->GetN() - 1, Fitting_ParamsError_MAP[type][energy][face][0], f2->GetParError(2));
                        // G_Calibration_MAP[energy]->AddPoint(Fitting_Params_MAP[type][energy][face][1], f2->GetParameter(4));
                        // G_Calibration_MAP[energy]->SetPointError(G_Calibration_MAP[energy]->GetN() - 1, Fitting_ParamsError_MAP[type][energy][face][1], f2->GetParError(4));
                        // G_Calibration_MAP[energy]->AddPoint(Fitting_Params_MAP[type][energy][face][2], f2->GetParameter(9));
                        // G_Calibration_MAP[energy]->SetPointError(G_Calibration_MAP[energy]->GetN() - 1, Fitting_ParamsError_MAP[type][energy][face][2], f2->GetParError(9));
                        // G_Calibration_MAP[energy]->AddPoint(Fitting_Params_MAP[type][energy][face][3], f2->GetParameter(11));
                        // G_Calibration_MAP[energy]->SetPointError(G_Calibration_MAP[energy]->GetN() - 1, Fitting_ParamsError_MAP[type][energy][face][5], f2->GetParError(11));
                        G_Calibration_MAP[energy]->AddPoint(Fitting_Params_MAP[type][energy][face][4], f2->GetParameter(14));
                        // G_Calibration_MAP[energy]->SetPointError(G_Calibration_MAP[energy]->GetN() - 1, Fitting_ParamsError_MAP[type][energy][face][6], f2->GetParError(14));
                    }
                }
            }
        }
    }

    for (const auto& energy : Energy)
    {
        if (energy == 3.0)
            continue;
        
        f->cd();
        G_Calibration_MAP[energy]->SetName(("G_Calibration_" + to_string(energy)).c_str());
        /// cout all points 
        for (int i = 0; i < G_Calibration_MAP[energy]->GetN(); i++)
        {
            double x, y;
            G_Calibration_MAP[energy]->GetPoint(i, x, y);
        }
        G_Calibration_MAP[energy]->Fit("pol1");
        TCanvas *c = new TCanvas("c", "c", 800, 800);
        G_Calibration_MAP[energy]->SetMarkerStyle(20);

        G_Calibration_MAP[energy]->Draw("AP");
        c->Write();

        cout << "Calibration Coefficients: " << endl;
        cout << "Offset 1.2 MeV: " << G_Calibration_MAP[energy]->GetFunction("pol1")->GetParameter(0) << endl;
        cout << "Coefficient 1.2 MeV: " << G_Calibration_MAP[energy]->GetFunction("pol1")->GetParameter(1) << endl;
    }
    

    

    for (const auto &type : Type)
    {
        // if (type == "THIN")
        //     continue;
        for (const auto &energy : Energy)
        {
            if (energy == 3.0)
                continue;
            for (const auto &face : Face)
            {
                int integral = Exp_Hist_MAP[type][energy][face]->Integral();

                //////////////////// ############### EXPERIMENTAL ############### ////////////////////
                Exp_Hist_calib = (TH1D *)Sim_Hist_conv_MAP[type][energy][face]->Clone("Exp_Hist_calib");
                Exp_Hist_calib->Reset();

                for (int bin = 0; bin < integral; bin++)
                {
                    Exp_Hist_calib->Fill(G_Calibration_MAP[energy]->GetFunction("pol1")->GetParameter(0) + G_Calibration_MAP[energy]->GetFunction("pol1")->GetParameter(1) * Exp_Hist_MAP[type][energy][face]->GetRandom());
                }
                Exp_Hist_calib_MAP[type][energy][face] = (TH1D *)Exp_Hist_calib->Clone(("Exp_Hist_calib_" + face).c_str());
                //////////////////// ############### chi2 ############### ////////////////////

                Exp_Hist_calib_MAP[type][energy][face]->GetXaxis()->SetRangeUser(fMIN_MAP[type][energy], fMAX_MAP[type][energy]);
                Sim_Hist_conv_MAP[type][energy][face]->GetXaxis()->SetRangeUser(fMIN_MAP[type][energy], fMAX_MAP[type][energy]);

                CHI2_MAP[type][energy][face] = Exp_Hist_calib_MAP[type][energy][face]->Chi2Test(Sim_Hist_conv_MAP[type][energy][face], "CHI2/NDF");
                // CHI2_MAP[type][energy][face] = Chi2Test(Exp_Hist_calib_MAP[type][energy][face], Sim_Hist_conv_MAP[type][energy][face], fMIN_MAP[type][energy], fMAX_MAP[type][energy]) / (fMAX_MAP[type][energy] - fMIN_MAP[type][energy]);

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

    CreateRootHist();
    FitExp();
    // Fitting_Params_MAP["THIN"][1.2]["A"] = {425.358, 445.156, 458.565, 479.365, 502.229, 525.009};
    // Fitting_ParamsError_MAP["THIN"][1.2]["A"] = {0.129836, 0.127157, 0.225123, 0.232471, 0.0511049, 0.0445305};
    // Fitting_Params_MAP["THIN"][1.2]["B"] = {418.947, 436.367, 453.849, 471.802, 499.497, 520.363};
    // Fitting_ParamsError_MAP["THIN"][1.2]["B"] = {0.136001, 0.134236, 0.251233, 0.24343, 0.0680359, 0.0474097};

    // Fitting_Params_MAP["THICK"][1.2]["A"] = {228.804, 443.989, 450, 200, 523.758};
    // Fitting_ParamsError_MAP["THICK"][1.2]["A"] = {0.357234, 0.097625, 5.17209, 0.0716736, 0.037395};

    // Fitting_Params_MAP["THICK"][1.2]["B"] = {188.683, 437.123, 440.988, 279.462, 520.46};
    // Fitting_ParamsError_MAP["THICK"][1.2]["B"] = {0.405468, 0.2703, 1.46313, 0.435669, 0.0415767};

    // INITIAL GUESS //


    // OFFSET
    FaceOffSet["A"] = "-2";
    FaceOffSet["B"] = "5";  

    // THICKNESS
    Thickness_Al1_MAP["THIN"] = 85;
    Thickness_Mylar_MAP["THIN"] = 525;
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
    file << chi2;
    file.close();

    
    double chi2_sum = 0;
    if (THREAD == 0)
    {
        
        f->cd();
        for (const auto &type : Type)
        {
            // if (type == "THIN")
            //     continue;
            for (const auto &energy : Energy)
            {
                if (energy == 3.0)
                    continue;
                for (const auto &face : Face)
                {
                    Exp_Hist_calib_MAP[type][energy][face]->Write();
                    Sim_Hist_conv_MAP[type][energy][face]->Write();

                    TCanvas *c = new TCanvas("c", "c", 800, 800);
                    Exp_Hist_calib_MAP[type][energy][face]->Draw("HIST");
                    Sim_Hist_conv_MAP[type][energy][face]->SetLineColor(kRed);
                    Sim_Hist_conv_MAP[type][energy][face]->Draw("HIST SAME");
                    cout << "- " << type << "\t" << energy << "\t" << face << "\t" << CHI2_MAP[type][energy][face] << endl;
                    chi2_sum += CHI2_MAP[type][energy][face];
                    c->Write();

                }
            }
        }
        
    }
    cout << "chi2_sum: " << chi2_sum << endl;
    f->Close();

    
    fSave->Close();

    return 0;
}