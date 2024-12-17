using namespace std;
using namespace ROOT::Math;

#include <iostream>
#include "TFile.h"
#include <fstream>
#include <dirent.h>
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"


vector<string> Face = {"A", "B"};
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

double fMIN = -1111;
double fMAX = -1111;

double energy;
TFile *Sim_File;
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
    string filename_f = "../../../../../../../mnt/hgfs/shared-2/ROOT_files/" + filename;
    system(("cd ../build && example " + macro_filename + " " + filename_f).c_str());

    return filename_f;
}

void CreateRootHist()
{
    for (string type : Type)
    {
        for (double energy : Energy)
        {
            for (string face : Face)
            {
                Exp_Hist_MAP[type][energy][face] = new TH1D(("Exp_Hist_" + type + "_" + to_string(energy) + "_" + face).c_str(), "Experimental RBS", 1024, 0, 1024);
                
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
                    if (line.compare(0, 14, "[DATA15,1024 ]") == 0)
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
    
    // Calibration_Offset_MAP[1.2] = par[0];
    Calibration_Offset_MAP[1.2] = 0;
    Calibration_Coefficient_MAP[1.2] = par[1];
    // Calibration_Offset_MAP[3.0] = par[2];
    Calibration_Offset_MAP[3.0] = 0;
    Calibration_Coefficient_MAP[3.0] = par[3];

    // round value to nearest 5 nm
    int step1 = 5;
    int step2 = 5; 
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
    for (string type : Type)
    {
        for (double energy : Energy)
        {
            for (string face : Face)
            {
                Exp_Hist = Exp_Hist_MAP[type][energy][face];
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
                    string new_macro_filename = MacroModifier(Macro_File, type, energy, face);
                    LaunchG4(new_macro_filename);
                    Sim_File = new TFile(filename.c_str(), "READ");
                }

                Al1 = (TH1D *)Sim_File->Get("RBS_0_13");
                Oxygen = (TH1D *)Sim_File->Get("RBS_1_8");
                Carbon = (TH1D *)Sim_File->Get("RBS_1_6");
                Al2 = (TH1D *)Sim_File->Get("RBS_2_13");

                Exp_Hist_calib = (TH1D *)Al1->Clone("Exp_Hist_calib");
                Exp_Hist_calib->Reset();

                int integral = Exp_Hist->Integral();

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

                Sim_Hist->GetXaxis()->SetRangeUser(Sim_WinMIN_MAP[type][energy], Sim_WinMAX_MAP[type][energy]);
                Exp_Hist->GetXaxis()->SetRangeUser(Exp_WinMIN_MAP[type][energy], Exp_WinMAX_MAP[type][energy]);
                Calibration_Offset_MAP[energy] += Sim_Hist->GetMean() - Calibration_Coefficient_MAP[energy] * Exp_Hist->GetMean();
                Sim_Hist->GetXaxis()->SetRangeUser(-1111, -1111);
                Exp_Hist->GetXaxis()->SetRangeUser(-1111, -1111);

                Sim_File->Close();
            }
        }
    }

    for (string type : Type)
    {
        for (double energy : Energy)
        {
            for (string face : Face)
            {
                Exp_Hist = Exp_Hist_MAP[type][energy][face];
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

                Al1 = (TH1D *)Sim_File->Get("RBS_0_13");
                Oxygen = (TH1D *)Sim_File->Get("RBS_1_8");
                Carbon = (TH1D *)Sim_File->Get("RBS_1_6");
                Al2 = (TH1D *)Sim_File->Get("RBS_2_13");

                Exp_Hist_calib = (TH1D *)Al1->Clone("Exp_Hist_calib");
                Exp_Hist_calib->Reset();

                int integral = Exp_Hist->Integral();

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

                for (int bin = 0; bin < integral; bin++)
                {
                    double value = Exp_Hist->GetRandom();
                    double value_calib = Calibration_Offset_MAP[energy]/4 + Calibration_Coefficient_MAP[energy] * value;
                    Exp_Hist_calib->Fill(value_calib);
                }

                // convolution sim_hist * gaussian
                TF1 *Gaussian = new TF1("Gaussian", "gausn", -1000, 1000);
                Gaussian->SetParameters(1, 0, 4);
                Sim_Hist_conv = (TH1D *)Sim_Hist->Clone("Sim_Hist_conv");
                Sim_Hist_conv->Reset();

                for (int bin = 0; bin < Sim_Hist->GetNbinsX(); bin++)
                {
                    double res = 0;
                    for (int bin_c = 0; bin_c < Exp_Hist_calib->GetNbinsX(); bin_c++)
                    {
                        res += Gaussian->Eval(Sim_Hist->GetBinCenter(bin) - Exp_Hist_calib->GetBinCenter(bin_c)) * Sim_Hist->GetBinContent(bin_c);
                    }
                    Sim_Hist_conv->SetBinContent(bin, res);
                    Sim_Hist_conv->SetBinError(bin, sqrt(res));
                }



                Al1_MAP[type][energy][face] = (TH1D *)Al1->Clone(("Al1_" + face).c_str());
                Oxygen_MAP[type][energy][face] = (TH1D *)Oxygen->Clone(("Oxygen_" + face).c_str());
                Carbon_MAP[type][energy][face] = (TH1D *)Carbon->Clone(("Carbon_" + face).c_str());
                Al2_MAP[type][energy][face] = (TH1D *)Al2->Clone(("Al2_" + face).c_str());
                Sim_Hist_MAP[type][energy][face] = (TH1D *)Sim_Hist->Clone(("Sim_Hist_" + face).c_str());
                Sim_Hist_conv_MAP[type][energy][face] = (TH1D *)Sim_Hist_conv->Clone(("Sim_Hist_conv_" + face).c_str());
                Exp_Hist_calib_MAP[type][energy][face] = (TH1D *)Exp_Hist_calib->Clone(("Exp_Hist_calib_" + face).c_str());

                Exp_Hist_calib_MAP[type][energy][face]->GetXaxis()->SetRangeUser(fMIN_MAP[type][energy], fMAX_MAP[type][energy]);
                Sim_Hist_conv_MAP[type][energy][face]->GetXaxis()->SetRangeUser(fMIN_MAP[type][energy], fMAX_MAP[type][energy]);

                CHI2_MAP[type][energy][face] = Exp_Hist_calib_MAP[type][energy][face]->Chi2Test(Sim_Hist_conv_MAP[type][energy][face], "CHI2/NDF");
                CHI2_SUM += CHI2_MAP[type][energy][face];

                Sim_File->Close();
            }
        }
    }

    // PRINTING CHI2
    cout << "CHI2: " << CHI2_SUM/8 << endl;
    cout << "Calibration Offset [THIN]: " << Calibration_Offset_MAP[1.2] << " keV" << endl;
    cout << "Calibration Coefficient [THIN]: " << Calibration_Coefficient_MAP[1.2] << " keV/CH" << endl;
    cout << "Calibration Offset [THICK]: " << Calibration_Offset_MAP[3.0] << " keV" << endl;
    cout << "Calibration Coefficient [THICK]: " << Calibration_Coefficient_MAP[3.0] << " keV/CH" << endl;

    cout << "Al1 Thickness [THIN]: " << Thickness_Al1_MAP["THIN"] << " nm" << endl;
    cout << "Mylar Thickness [THIN]: " << Thickness_Mylar_MAP["THIN"] << " nm" << endl;
    cout << "Al2 Thickness [THIN]: " << Thickness_Al2_MAP["THIN"] << " nm" << endl;

    cout << "Al1 Thickness [THICK]: " << Thickness_Al1_MAP["THICK"] << " nm" << endl;
    cout << "Mylar Thickness [THICK]: " << Thickness_Mylar_MAP["THICK"] << " nm" << endl;
    cout << "Al2 Thickness [THICK]: " << Thickness_Al2_MAP["THICK"] << " nm" << endl;
    for (string type : Type)
    {
        for (double energy : Energy)
        {
            for (string face : Face)
            {
                if (CHI2_MAP[type][energy][face] > 10)
                {
                    cout << "\033[1;31m";
                }
                cout << "- " << type << "\t" << energy << "\t" << face << "\t" << CHI2_MAP[type][energy][face] << "\033[0m" << endl;
            }
        }
    }
    cout << endl;

    return CHI2_SUM;
}

void RBS_All()
{

    TFile *f = new TFile("RBS_Results.root", "RECREATE");


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

    // CALIBRATION
    // Calibration_Offset_MAP[1.2] = 1.89;
    // Calibration_Coefficient_MAP[1.2] = 2.00495;

    // Calibration_Offset_MAP[3.0] = 25;
    // Calibration_Coefficient_MAP[3.0] = 3.7989;

    //moi
    Calibration_Offset_MAP[1.2] = 85.2619;
    Calibration_Coefficient_MAP[1.2] = 1.9697;

    Calibration_Offset_MAP[3.0] = 288.154;
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
    


    Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Functor functor(&FunctionToMinimize, 10);
    minimizer->SetFunction(functor);
    // CALIBRATION PARAMETERS
    // minimizer->SetFixedVariable(0, "Calibration Offset [THIN]", Calibration_Offset_MAP[1.2]);
    minimizer->SetLimitedVariable(0, "Calibration Offset [THIN]", Calibration_Offset_MAP[1.2], 20, 0, 200);

    // minimizer->SetFixedVariable(1, "Calibration Coefficient [THIN]", Calibration_Coefficient_MAP[1.2]);
    minimizer->SetLimitedVariable(1, "Calibration Coefficient [THIN]", Calibration_Coefficient_MAP[1.2], 0.1, 0.9 * Calibration_Coefficient_MAP[1.2], 2.01);

    // minimizer->SetFixedVariable(2, "Calibration Offset [THICK]", Calibration_Offset_MAP[3.0]);
    minimizer->SetLimitedVariable(2, "Calibration Offset [THICK]", Calibration_Offset_MAP[3.0], 20, 0, 500);

    //  minimizer->SetFixedVariable(3, "Calibration Coefficient [THICK]", Calibration_Coefficient_MAP[3.0]);
    minimizer->SetLimitedVariable(3, "Calibration Coefficient [THICK]", Calibration_Coefficient_MAP[3.0], 0.2, 0.9 * Calibration_Coefficient_MAP[3.0], 3.8);  


    // THICKNESS PARAMETERS
    // minimizer->SetFixedVariable(4, "Al1 Thickness [THIN]", Thickness_Al1_MAP["THIN"]);
    minimizer->SetLimitedVariable(4, "Al1 Thickness [THIN]", Thickness_Al1_MAP["THIN"], 20, 50, 200);
    // minimizer->SetFixedVariable(5, "Mylar Thickness [THIN]", Thickness_Mylar_MAP["THIN"]);
    minimizer->SetLimitedVariable(5, "Mylar Thickness [THIN]", Thickness_Mylar_MAP["THIN"], 10, 400, 700);
    minimizer->SetFixedVariable(6, "Al2 Thickness [THIN]", Thickness_Al2_MAP["THIN"]);
    // minimizer->SetLimitedVariable(6, "Al2 Thickness [THIN]", Thickness_Al2_MAP["THIN"], 20, 50, 200);
    // minimizer->SetFixedVariable(7, "Al1 Thickness [THICK]", Thickness_Al1_MAP["THICK"]);
    minimizer->SetLimitedVariable(7, "Al1 Thickness [THICK]", Thickness_Al1_MAP["THICK"], 20, 50, 200);
    // minimizer->SetFixedVariable(8, "Mylar Thickness [THICK]", Thickness_Mylar_MAP["THICK"]);
    minimizer->SetLimitedVariable(8, "Mylar Thickness [THICK]", Thickness_Mylar_MAP["THICK"], 10, 5500, 6500);
    minimizer->SetFixedVariable(9, "Al2 Thickness [THICK]", Thickness_Al2_MAP["THICK"]);
    // minimizer->SetLimitedVariable(9, "Al2 Thickness [THICK]", Thickness_Al2_MAP["THICK"], 20, 50, 200);

    // SCALE PARAMETERS
    int counter = 10;
    // for (string type : Type)
    // {
    //     for (double energy : Energy)
    //     {
    //         minimizer->SetFixedVariable(counter, ("Al1 Scale " + type + " " + to_string(energy)).c_str(), ScaleAl1_MAP[type][energy]);
    //         minimizer->SetFixedVariable(counter + 1, ("Oxygen Scale " + type + " " + to_string(energy)).c_str(), ScaleOxygen_MAP[type][energy]);
    //         minimizer->SetFixedVariable(counter + 2, ("Carbon Scale " + type + " " + to_string(energy)).c_str(), ScaleCarbon_MAP[type][energy]);
    //         minimizer->SetFixedVariable(counter + 3, ("Al2 Scale " + type + " " + to_string(energy)).c_str(), ScaleAl2_MAP[type][energy]);
    //         counter += 4;
    //     }
    // }
    
    // minimizer->SetPrecision(0.001);
    // minimizer->SetTolerance(0.001);
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetMaxIterations(1000000);
    minimizer->SetPrintLevel(1);

    minimizer->Minimize();
    const double *par = minimizer->X();

    // const double par[10] = {
    //     Calibration_Offset_MAP[1.2], Calibration_Coefficient_MAP[1.2],
    //     Calibration_Offset_MAP[3.0], Calibration_Coefficient_MAP[3.0],
    //     Thickness_Al1_MAP["THIN"], Thickness_Mylar_MAP["THIN"], Thickness_Al2_MAP["THIN"],
    //     Thickness_Al1_MAP["THICK"], Thickness_Mylar_MAP["THICK"], Thickness_Al2_MAP["THICK"]
    // };

    //// GETCOVARAINCE MATRIX
    // TMatrixDSym cov = minimizer->CovMatrix();
    // cout << "Covariance Matrix: " << endl;
    // cov.Print();

    //// ERRORS from method Errors()
    const double *err = minimizer->Errors();
    cout << "Errors: " << endl;
    for (int i = 0; i < 10; i++)
    {
        cout << err[i] << endl;
    }


    FunctionToMinimize(par);

    // gStyle->SetOptStat(0);
    
    f->cd();

    TCanvas *c = new TCanvas("All", "All", 800, 600);
    c->Divide(2, 4);
    counter = 1;
    for (string type : Type)
    {
        for (double energy : Energy)
        {
            ostringstream oss;
            oss << fixed << setprecision(1) << energy;

            TCanvas *c1 = new TCanvas(("CHI2_" + type + "_" + oss.str()).c_str(), ("CHI2_" + type + "_" + oss.str()).c_str(), 800, 600);
            c1->Divide(1, 2);
            TLegend *legend_1 = new TLegend(0.7, 0.7, 0.9, 0.9);
            TCanvas *c2 = new TCanvas(("Components_" + type + "_" + oss.str()).c_str(), ("Components_" + type + "_" + oss.str()).c_str(), 800, 600);
            c2->Divide(2, 1);
            TLegend *legend_2 = new TLegend(0.7, 0.7, 0.9, 0.9);
            int face_counter = 1;
            for (string face : Face)
            {
                gStyle->SetOptStat(0);
                ostringstream osss;
                osss << fixed << setprecision(2) << CHI2_MAP[type][energy][face];
                TLatex *t = new TLatex(fMIN_MAP[type][energy] + (fMAX_MAP[type][energy] - fMIN_MAP[type][energy]) * 0.8, Exp_Hist_calib_MAP[type][energy][face]->GetMaximum() * 0.65, ("#chi^{2}_{#nu} = " + osss.str()).c_str());
                t->SetTextSize(0.1);
                t->SetTextColor(kRed);

                c->cd(counter);
                Exp_Hist_calib_MAP[type][energy][face]->SetLineColor(kBlack);
                Exp_Hist_calib_MAP[type][energy][face]->GetXaxis()->SetTitle("Energy [keV]");
                Exp_Hist_calib_MAP[type][energy][face]->GetYaxis()->SetTitle("Counts/keV");
                Exp_Hist_calib_MAP[type][energy][face]->SetTitle((type + " " + oss.str() + " MeV " + face).c_str());
                Exp_Hist_calib_MAP[type][energy][face]->Draw("HIST");

                Sim_Hist_conv_MAP[type][energy][face]->SetLineColor(kRed);
                Sim_Hist_conv_MAP[type][energy][face]->Draw("HIST SAME");

                t->Draw("SAME");

                if (counter == 1)
                {
                    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
                    legend->AddEntry(Exp_Hist_calib_MAP[type][energy][face], "Experimental", "l");
                    legend->AddEntry(Sim_Hist_conv_MAP[type][energy][face], "Simulated", "l");
                    legend->Draw("SAME");
                }
                counter++;


                c1->cd(face_counter);
                Exp_Hist_calib_MAP[type][energy][face]->Draw("HIST");
                legend_1->AddEntry(Exp_Hist_calib_MAP[type][energy][face], "Experimental", "l");
                Sim_Hist_conv_MAP[type][energy][face]->Draw("HIST SAME");
                legend_1->AddEntry(Sim_Hist_conv_MAP[type][energy][face], "Simulated", "l");
                t->Draw("SAME");
                                
                c2->cd(face_counter);
                Exp_Hist_calib_MAP[type][energy][face]->Draw("HIST");
                legend_2->AddEntry(Exp_Hist_calib_MAP[type][energy][face], "Experimental", "l");
                Al1_MAP[type][energy][face]->SetLineColor(kCyan);
                Al1_MAP[type][energy][face]->SetTitle("Al (front)");
                Al1_MAP[type][energy][face]->Draw("HIST SAME");
                legend_2->AddEntry(Al1_MAP[type][energy][face], "Al (front)", "l");
                Al2_MAP[type][energy][face]->SetLineColor(kGreen);
                Al2_MAP[type][energy][face]->SetTitle("Al (back)");
                Al2_MAP[type][energy][face]->Draw("HIST SAME");
                legend_2->AddEntry(Al2_MAP[type][energy][face], "Al (back)", "l");
                Oxygen_MAP[type][energy][face]->SetLineColor(kBlue);
                Oxygen_MAP[type][energy][face]->SetTitle("Oxygen");
                Oxygen_MAP[type][energy][face]->Draw("HIST SAME");
                legend_2->AddEntry(Oxygen_MAP[type][energy][face], "Oxygen", "l");
                Carbon_MAP[type][energy][face]->SetLineColor(kMagenta);
                Carbon_MAP[type][energy][face]->SetTitle("Carbon");
                Carbon_MAP[type][energy][face]->Draw("HIST SAME");
                legend_2->AddEntry(Carbon_MAP[type][energy][face], "Carbon", "l");
                face_counter++;

                Exp_Hist_MAP[type][energy][face]->SetTitle((type + " " + oss.str() + " MeV " + face).c_str());
                Exp_Hist_MAP[type][energy][face]->Write();

            }
            c1->cd();
            legend_1->Draw("SAME");
            c1->Write();
            c2->cd();
            legend_2->Draw("SAME");
            c2->Write();

            c1->Close();
            c2->Close();
        }
    }
    c->Draw();
    c->Write();

    f->Close();
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

