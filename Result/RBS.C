using namespace std;
using namespace  ROOT::Math;

#include <iostream>
#include "TFile.h"
#include <fstream>
#include <dirent.h>
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

TH1D *Exp_Hist = new TH1D("Exp_Hist", "Experimental RBS", 1024, 0, 1024);
TH1D *Al1;
TH1D *Oxygen;
TH1D *Carbon;
TH1D *Al2;

TH1D *Exp_Hist_calib;
TH1D *Sim_Hist;
TH1D* Sim_Hist_conv;

double energy;
TFile *Sim_File;
string Macro_File;

string MacroModifier(string filename, const double *par, double energy)
{
    double T1 = par[2];
    double T2 = par[3];
    double T3 = par[4];

    ostringstream oss;
    oss << fixed << setprecision(1) << energy;
    string new_filename = "../shoot_" + oss.str() + "_" + to_string((int)T1) + "_" + to_string((int)T2) + "_" + to_string((int)T3) + ".mac";

    ifstream file;
    file.open("../shoot_base.mac");
    ofstream temp_file(new_filename);
    string line;

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
        temp_file << line << endl;
    }

    file.close();
    temp_file.close();

    return new_filename;
}

string LaunchG4(string macro_filename)
{
    // 9 to end of string
    string filename = macro_filename.substr(9, macro_filename.size() - 13) + ".root";
    system(("cd ../build && example " + macro_filename + " " + filename).c_str());

    return filename;
}

void CreateRootHist(string filename)
{
    // read txt file
    ifstream file;
    file.open(filename);
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
            Exp_Hist->SetBinContent(counter, value);
            counter++;
        }
        if (line.compare(0, 1, "[") == 0)
        {
            data = false;
        }
    }

    file.close();
}

double FunctionToMinimize(const double *par)
{

    // loop
    // minimizer

    double b = par[0];
    double a = par[1];
    double thickness_Al1 = par[2];
    double thickness_Mylar = par[3];
    double thickness_Al2 = par[4];
    double scaleA1 = par[5];
    double scaleOxygen = par[6];
    double scaleCarbon = par[7];
    double scaleA2 = par[8];

    ostringstream oss;
    oss << fixed << setprecision(1) << energy;
    string ss = oss.str();
    string filename = oss.str() + "_" + to_string((int)thickness_Al1) + "_" + to_string((int)thickness_Mylar) + "_" + to_string((int)thickness_Al2) + ".root";
    Sim_File = new TFile(filename.c_str(), "READ");
    if (Sim_File->IsZombie())
    {
        string new_macro_filename = MacroModifier(Macro_File, par, energy);
        string filename = LaunchG4(new_macro_filename);
        Sim_File = new TFile(filename.c_str(), "READ");
    }

    Al1 = (TH1D *)Sim_File->Get("RBS_0_13");
    Oxygen = (TH1D *)Sim_File->Get("RBS_1_8");
    Carbon = (TH1D *)Sim_File->Get("RBS_1_6");
    Al2 = (TH1D *)Sim_File->Get("RBS_2_13");

    Exp_Hist_calib = (TH1D *)Al1->Clone("Exp_Hist_calib");
    Exp_Hist_calib->Reset();

    int integral = Exp_Hist->Integral();

    Exp_Hist->GetXaxis()->SetRangeUser(650, 750);
    Al1->GetXaxis()->SetRangeUser(2500, 2700);
    b = 18.25;//Al1->GetMean()-a*Exp_Hist->GetMean();

    Exp_Hist->GetXaxis()->SetRangeUser(-1111, -1111);
    Al1->GetXaxis()->SetRangeUser(-1111, -1111);


    for (int bin = 0; bin < integral; bin++)
    {
        double value = Exp_Hist->GetRandom();
        double value_calib = b + a * value;
        Exp_Hist_calib->Fill(value_calib);
    }

    Al1->Scale(scaleA1 / Al1->Integral());
    Oxygen->Scale(scaleOxygen / Oxygen->Integral());
    Carbon->Scale(scaleCarbon / Carbon->Integral());
    Al2->Scale(scaleA2 / Al2->Integral());

    Sim_Hist = (TH1D *)Al1->Clone("Sim_Hist");
    Sim_Hist->Reset();
    Sim_Hist->Add(Al1);
    Sim_Hist->Add(Oxygen);
    Sim_Hist->Add(Carbon);
    Sim_Hist->Add(Al2);

    // convolution sim_hist * gaussian
    TF1 *Gaussian = new TF1("Gaussian", "gausn", -1000, 1000);
    Gaussian->SetParameters(1, 0, 3.5);
    Sim_Hist_conv = (TH1D*)Sim_Hist->Clone("Sim_Hist_conv");
    Sim_Hist_conv->Reset();

    for (int bin = 0; bin < Sim_Hist->GetNbinsX(); bin++)
    {
        double res = 0;
        for (int bin_c = 0; bin_c < Exp_Hist_calib->GetNbinsX(); bin_c++)
        {
            res += Gaussian->Eval(Sim_Hist->GetBinCenter(bin) - Exp_Hist_calib->GetBinCenter(bin_c)) * Sim_Hist->GetBinContent(bin_c);
        }
        Sim_Hist_conv->SetBinContent(bin, res);
    }


    double chi2 = Exp_Hist_calib->Chi2Test(Sim_Hist, "CHI2/NDF");

    cout << chi2 << endl;

    return chi2;
}

void RBS()
{

    // random generator
    random_device rd;
    mt19937 generator(rd());

    // string Exp_FileName = "2023W43/RBS_Pos1.mpa";
    string Exp_FileName = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_A_Pos2_3MeV.mpa";
    energy = 3.0; // MeV

    
    Macro_File = "../shoot.mac";

    CreateRootHist(Exp_FileName);

    double thickness_Al1 = 116;   // nm
    double thickness_Mylar = 6200; // nm
    double thickness_Al2 = 106;   // nm

    double scaleA1 = 1000;
    double scaleOxygen = 10000;
    double scaleCarbon = 10000;
    double scaleA2 = 1000;

    double a = 1;
    if (energy == 1.2)
    {
        a = 2;
    }
    else
    {
        a = 3.6;
    }

    Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Functor functor(&FunctionToMinimize, 9);
    minimizer->SetFunction(functor);
    minimizer->SetFixedVariable(0, "Calibration Offset", 0);
    minimizer->SetLimitedVariable(1, "Calibration Coefficient", a, 0.01, 0.9 * a, 1.1 * a);
    minimizer->SetFixedVariable(2, "Al1 Thickness", thickness_Al1);
    minimizer->SetFixedVariable(3, "Mylar Thickness", thickness_Mylar);
    minimizer->SetFixedVariable(4, "Al2 Thickness", thickness_Al2);
    //1.2
    // minimizer->SetFixedVariable(5, "Al1 Scale", scaleA1*5);
    // minimizer->SetFixedVariable(6, "Oxygen Scale", scaleOxygen*6);
    // minimizer->SetFixedVariable(7, "Carbon Scale", scaleCarbon*22);
    // minimizer->SetFixedVariable(8, "Al2 Scale", scaleA2*13);
    //3.0
    minimizer->SetFixedVariable(5, "Al1 Scale", scaleA1*3.8);
    minimizer->SetFixedVariable(6, "Oxygen Scale", scaleOxygen*6.4);
    minimizer->SetFixedVariable(7, "Carbon Scale", scaleCarbon*16.2);
    minimizer->SetFixedVariable(8, "Al2 Scale", scaleA2*3.5);
    minimizer->SetPrecision(0.001);
    minimizer->SetTolerance(0.001);
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetMaxIterations(1000000);

    // minimizer->Minimize();
    // const double *par = minimizer->X();

    // const double par[9] = {-20, 2.05, thickness_Al1, thickness_Mylar, thickness_Al2, scaleA1*5, scaleOxygen*6, scaleCarbon*22, scaleA2*13};

    // const double par[9] = {-28.76, 2.05, thickness_Al1, thickness_Mylar, thickness_Al2, scaleA1*5, scaleOxygen*6, scaleCarbon*22, scaleA2*13};



    //RBS_A_Pos2_3MeV
    //116 6400 106
    // const double par[9] = {0, 3.79147, thickness_Al1, thickness_Mylar, thickness_Al2, scaleA1*3.7, scaleOxygen*6.4, scaleCarbon*16.2, scaleA2*1.5};
    const double par[9] = {0, 3.81147, thickness_Al1, thickness_Mylar, thickness_Al2, scaleA1*3.7, scaleOxygen*6.4, scaleCarbon*16.2, scaleA2*1.5};

    //RBS_B_Pos2_3MeV
    //116 6400 106
    double coef= 121546./238556.;
    // const double par[9] = {0, 3.79147, thickness_Al1, thickness_Mylar, thickness_Al2, scaleA1*3.7*coef, scaleOxygen*6.4*coef, scaleCarbon*16.2*coef, scaleA2*1.5*coef};

    // print par 
    cout << "Calibration Offset: " << par[0] << endl;
    cout << "Calibration Coefficient: " << par[1] << endl;


    FunctionToMinimize(par);

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    Exp_Hist_calib->SetLineColor(kBlack);
    Exp_Hist_calib->SetTitle("Experimental");
    Exp_Hist_calib->Draw("HIST");
    // Al1->SetLineColor(kCyan);
    // Al1->SetTitle("Al (front)");
    // Al1->Draw("HIST SAME");
    // Oxygen->SetLineColor(kBlue);
    // Oxygen->SetTitle("Oxygen");
    // Oxygen->Draw("HIST SAME");
    // Carbon->SetLineColor(kMagenta);
    // Carbon->SetTitle("Carbon");
    // Carbon->Draw("HIST SAME");
    // Al2->SetLineColor(kGreen);
    // Al2->SetTitle("Al (back)");
    // Al2->Draw("HIST SAME");
    // Sim_Hist->SetLineColor(kRed);
    // Sim_Hist->SetTitle("Simulated");
    // Sim_Hist->Draw("HIST SAME");
    Sim_Hist_conv->SetLineColor(kOrange);
    Sim_Hist_conv->SetLineWidth(2);
    Sim_Hist_conv->SetTitle("Simulated Conv");
    Sim_Hist_conv->Draw("HIST SAME");
    c->BuildLegend();
    c->Draw();
}