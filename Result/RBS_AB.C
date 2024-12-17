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
vector<int> FaceOffSet = {-2, 6};

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

string MacroModifier(string filename, const double *par, double energy, string Face)
{
    double T1 = par[2];
    double T2 = par[3];
    double T3 = par[4];

    ostringstream oss;
    oss << fixed << setprecision(1) << energy;
    string new_filename = "../tmp/shoot_" + oss.str() + "_" + to_string((int)T1) + "_" + to_string((int)T2) + "_" + to_string((int)T3) + "_" + Face + ".mac";

    ifstream file;
    file.open("../shoot_base.mac");
    ofstream temp_file(new_filename);
    string line;

    int face_index;
    if (Face == "A")
    {
        face_index = 0;
    }
    else
    {
        face_index = 1;
        // SWAP Al thickness
        T1 = par[4];
        T3 = par[2];
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
            line.replace(line.find("%p"), 2, to_string((int)FaceOffSet[face_index]) + " mm");
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
    string filename = macro_filename.substr(13, macro_filename.size() - 17) + ".root";
    system(("cd ../build && example " + macro_filename + " " + filename).c_str());

    return filename;
}

void CreateRootHist(vector<string> filename_pair)
{

    ////// A
    string filename = filename_pair[0];
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

    Exp_Hist_AB.push_back((TH1D *)Exp_Hist->Clone("Exp_Hist_A"));

    ////// B
    Exp_Hist->Reset();
    filename = filename_pair[1];
    // read txt file
    file.open(filename);
    data = false;
    counter = 0;

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

    Exp_Hist_AB.push_back((TH1D *)Exp_Hist->Clone("Exp_Hist_B"));
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

    double coef = 1;
    
    
    for (int i = 0; i < 2; i++)
    {
        Exp_Hist = Exp_Hist_AB[i];
        if (i == 1) // coef from experimental normalization
        {
            coef = (double)Exp_Hist_AB[1]->Integral() / (double)Exp_Hist_AB[0]->Integral();
        }

        string filename = oss.str() + "_" + to_string((int)thickness_Al1) + "_" + to_string((int)thickness_Mylar) + "_" + to_string((int)thickness_Al2) + "_" + Face[i] + ".root";
        Sim_File = new TFile(filename.c_str(), "READ");
        if (Sim_File->IsZombie())
        {
            string new_macro_filename = MacroModifier(Macro_File, par, energy, Face[i]);
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

        Exp_Hist->GetXaxis()->SetRangeUser(-1111, -1111);
        Al1->GetXaxis()->SetRangeUser(-1111, -1111);

        for (int bin = 0; bin < integral; bin++)
        {
            double value = Exp_Hist->GetRandom();
            double value_calib = b + a * value;
            Exp_Hist_calib->Fill(value_calib);
        }

        Al1->Scale(scaleA1*coef / Al1->Integral());
        Oxygen->Scale(scaleOxygen*coef / Oxygen->Integral());
        Carbon->Scale(scaleCarbon*coef / Carbon->Integral());
        Al2->Scale(scaleA2*coef / Al2->Integral());

        Sim_Hist = (TH1D *)Al1->Clone("Sim_Hist");
        Sim_Hist->Reset();
        Sim_Hist->Add(Al1);
        Sim_Hist->Add(Oxygen);
        Sim_Hist->Add(Carbon);
        Sim_Hist->Add(Al2);

        

        // convolution sim_hist * gaussian
        TF1 *Gaussian = new TF1("Gaussian", "gausn", -1000, 1000);
        Gaussian->SetParameters(1, 0, 3.5);
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

        Al1_AB.push_back((TH1D *)Al1->Clone(("Al1_" + Face[i]).c_str()));
        Oxygen_AB.push_back((TH1D *)Oxygen->Clone(("Oxygen_" + Face[i]).c_str()));
        Carbon_AB.push_back((TH1D *)Carbon->Clone(("Carbon_" + Face[i]).c_str()));
        Al2_AB.push_back((TH1D *)Al2->Clone(("Al2_" + Face[i]).c_str()));
        Sim_Hist_AB.push_back((TH1D *)Sim_Hist->Clone(("Sim_Hist_" + Face[i]).c_str()));
        Sim_Hist_conv_AB.push_back((TH1D *)Sim_Hist_conv->Clone(("Sim_Hist_conv_" + Face[i]).c_str()));
        Exp_Hist_calib_AB.push_back((TH1D *)Exp_Hist_calib->Clone(("Exp_Hist_calib_" + Face[i]).c_str()));
    }

    Exp_Hist_calib_AB[0]->GetXaxis()->SetRangeUser(fMIN, fMAX);
    Exp_Hist_calib_AB[1]->GetXaxis()->SetRangeUser(fMIN, fMAX);
    Sim_Hist_conv_AB[0]->GetXaxis()->SetRangeUser(fMIN, fMAX);
    Sim_Hist_conv_AB[1]->GetXaxis()->SetRangeUser(fMIN, fMAX);
    double chi2_0 = Exp_Hist_calib_AB[0]->Chi2Test(Sim_Hist_conv_AB[0], "CHI2/NDF");
    double chi2_1 = Exp_Hist_calib_AB[1]->Chi2Test(Sim_Hist_conv_AB[1], "CHI2/NDF");

    cout << a << "    " << b << endl;
    cout << (chi2_0 + chi2_1)/2 << "   (Chi2 A: " << chi2_0 << "      Chi2 B: " << chi2_1 << ")" << endl;

    //cleaning
    Al1_AB.clear();
    Oxygen_AB.clear();
    Carbon_AB.clear();
    Al2_AB.clear();
    Sim_Hist_AB.clear();
    Sim_Hist_conv_AB.clear();
    Exp_Hist_calib_AB.clear();


    
    return (chi2_0 + chi2_1)/2;
}

void RBS_AB()
{
    energy = 3.0; // MeV
    string type = "thin";
    vector<string> Exp_FileName = vector<string>(2);
    Exp_FileName[0] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_A_Pos1_1,2MeV.mpa";
    Exp_FileName[1] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_B_Pos1_1,2MeV.mpa";
    if (type == "thick")
    {
        if (energy == 3.0)
        {
            Exp_FileName[0] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_A_Pos2_3MeV.mpa";
            Exp_FileName[1] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_B_Pos2_3MeV.mpa";
            fMIN = 1900;
            fMAX = 3000;
        }
    
        if (energy == 1.2)
        {
            Exp_FileName[0] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_A_Pos2_1,2MeV.mpa";
            Exp_FileName[1] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_B_Pos2_1,2MeV.mpa";
            fMIN = 200;
            fMAX = 1200;
        }
    }
    else
    {
        if (energy == 3.0)
        {
            Exp_FileName[0] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_A_Pos1_3MeV.mpa";
            Exp_FileName[1] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_B_Pos1_3MeV.mpa";
            fMIN = 2100;
            fMAX = 2800;
        }

        if (energy == 1.2)
        {
            Exp_FileName[0] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_A_Pos1_1,2MeV.mpa";
            Exp_FileName[1] = "../../../../../../mnt/hgfs/shared-2/2024W35/RBS_B_Pos1_1,2MeV.mpa";
            fMIN = 800;
            fMAX = 1200;
        }
    }
    

    Macro_File = "../shoot.mac";

    CreateRootHist(Exp_FileName);

    double thickness_Al1 = 67;    // nm
    double thickness_Mylar = 600; // nm
    double thickness_Al2 = 50;    // nm

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
        a = 3.79147;
    }

    Minimizer *minimizer = Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Functor functor(&FunctionToMinimize, 9);
    minimizer->SetFunction(functor);
    // minimizer->SetVariableLimits(0, "Calibration Offset", 50, 0, 200);
    // minimizer->SetVariableLimits(1, "Calibration Coefficient", a, 0.9 * a, 1.1 * a);
    // minimizer->SetVariableLimits(2, "Al1 Thickness", thickness_Al1, 50, 200);
    // minimizer->SetVariableStepSize(2, 5);
    // minimizer->SetVariableLimits(3, "Mylar Thickness", thickness_Mylar, 400, 600);
    // minimizer->SetVariableStepSize(3, 10);
    // minimizer->SetVariableLimits(4, "Al2 Thickness", thickness_Al2, 50, 200);
    // minimizer->SetVariableStepSize(4, 5);

    minimizer->SetLimitedVariable(0, "Calibration Offset", 20, 20, 0, 200);
    minimizer->SetLimitedVariable(1, "Calibration Coefficient", a, 0.1, 0.9 * a, 1.1 * a);

    // minimizer->SetLimitedVariable(2, "Al1 Thickness", thickness_Al1, 10, 50, 200);
    // minimizer->SetLimitedVariable(3, "Mylar Thickness", thickness_Mylar, 20, 400, 600);
    // minimizer->SetLimitedVariable(4, "Al2 Thickness", thickness_Al2, 10, 50, 200);

    minimizer->SetFixedVariable(2, "Al1 Thickness", thickness_Al1);
    minimizer->SetFixedVariable(3, "Mylar Thickness", thickness_Mylar);
    minimizer->SetFixedVariable(4, "Al2 Thickness", thickness_Al2);


    if (type == "thin")
    {
        if (energy == 1.2)
        {
            // minimizer->SetLimitedVariable(5, "Al1 Scale", scaleA1 * 5, 100, scaleA1 * (5-0.2), scaleA1 * (5+0.2));
            // minimizer->SetLimitedVariable(6, "Oxygen Scale", scaleOxygen * 0.7, 100, scaleOxygen * 0.6, scaleOxygen * 0.8);
            // minimizer->SetLimitedVariable(7, "Carbon Scale", scaleCarbon * 2.3, 100, scaleCarbon * 2.2, scaleCarbon * 2.4);
            // minimizer->SetLimitedVariable(8, "Al2 Scale", scaleA2 * 5, 100, scaleA2 * 4.8, scaleA2 * 5.2);

            minimizer->SetFixedVariable(5, "Al1 Scale", scaleA1 * 5);
            minimizer->SetFixedVariable(6, "Oxygen Scale", scaleOxygen * 0.7);
            minimizer->SetFixedVariable(7, "Carbon Scale", scaleCarbon * 2.3);
            minimizer->SetFixedVariable(8, "Al2 Scale", scaleA2 * 5);
        }
        else
        {
            minimizer->SetFixedVariable(5, "Al1 Scale", scaleA1 * 2.1);
            minimizer->SetFixedVariable(6, "Oxygen Scale", scaleOxygen * 0.6);
            minimizer->SetFixedVariable(7, "Carbon Scale", scaleCarbon * 1.5);
            minimizer->SetFixedVariable(8, "Al2 Scale", scaleA2 * 1.2);
        }
    }
    else
    {
        if (energy == 1.2)
        {
            minimizer->SetFixedVariable(5, "Al1 Scale", scaleA1 * 8);
            minimizer->SetFixedVariable(6, "Oxygen Scale", scaleOxygen * 8);
            minimizer->SetFixedVariable(7, "Carbon Scale", scaleCarbon * 30);
            minimizer->SetFixedVariable(8, "Al2 Scale", scaleA2 * 13);
        }
        else
        {
            minimizer->SetFixedVariable(5, "Al1 Scale", scaleA1 * 3.7);
            minimizer->SetFixedVariable(6, "Oxygen Scale", scaleOxygen * 6.4);
            minimizer->SetFixedVariable(7, "Carbon Scale", scaleCarbon * 16.2);
            minimizer->SetFixedVariable(8, "Al2 Scale", scaleA2 * 1.5);
        }
    }
    minimizer->SetPrecision(0.001);
    minimizer->SetTolerance(0.001);
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetMaxIterations(1000000);

    minimizer->Minimize();
    const double *par = minimizer->X();

    // cout paramaters
    cout << "Calibration Offset: " << par[0] << endl;
    cout << "Calibration Coefficient: " << par[1] << endl;
    cout << "Al1 Thickness: " << par[2] << endl;
    cout << "Mylar Thickness: " << par[3] << endl;
    cout << "Al2 Thickness: " << par[4] << endl;
    cout << "Al1 Scale: " << par[5] << endl;
    cout << "Oxygen Scale: " << par[6] << endl;
    cout << "Carbon Scale: " << par[7] << endl;
    cout << "Al2 Scale: " << par[8] << endl;


    // const double par[9] = {-20, 2.05, thickness_Al1, thickness_Mylar, thickness_Al2, scaleA1*5, scaleOxygen*6, scaleCarbon*22, scaleA2*13};

    // const double par[9] = {-28.76, 2.05, thickness_Al1, thickness_Mylar, thickness_Al2, scaleA1*5, scaleOxygen*6, scaleCarbon*22, scaleA2*13};

    // RBS_A_Pos2_3MeV
    // 116 6400 106
    //  const double par[9] = {0, 3.79147, thickness_Al1, thickness_Mylar, thickness_Al2, scaleA1*3.7, scaleOxygen*6.4, scaleCarbon*16.2, scaleA2*1.5};
    // const double par[9] = {29.29, 3.79316, thickness_Al1, thickness_Mylar, thickness_Al2, scaleA1 * 2.1, scaleOxygen * 0.6, scaleCarbon * 1.5, scaleA2 * 1.2};
    
    // print par
    cout << "Calibration Offset: " << par[0] << endl;
    cout << "Calibration Coefficient: " << par[1] << endl;

    FunctionToMinimize(par);

    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    c->Divide(1, 2);
    for (int i = 0; i < 2; i++)
    {
        c->cd(i+1);
        Exp_Hist_calib_AB[i]->SetLineColor(kBlack);
        Exp_Hist_calib_AB[i]->SetTitle("Experimental");
        Exp_Hist_calib_AB[i]->Draw("HIST");

        // Al1_AB[i]->SetLineColor(kCyan);
        // Al1_AB[i]->SetTitle("Al (front)");
        // Al1_AB[i]->Draw("HIST SAME");
        // Oxygen_AB[i]->SetLineColor(kBlue);
        // Oxygen_AB[i]->SetTitle("Oxygen");
        // Oxygen_AB[i]->Draw("HIST SAME");
        // Carbon_AB[i]->SetLineColor(kMagenta);
        // Carbon_AB[i]->SetTitle("Carbon");
        // Carbon_AB[i]->Draw("HIST SAME");
        // Al2_AB[i]->SetLineColor(kGreen);
        // Al2_AB[i]->SetTitle("Al (back)");
        // Al2_AB[i]->Draw("HIST SAME");

        // Sim_Hist_AB[i]->SetLineColor(kOrange);
        // Sim_Hist_AB[i]->SetTitle("Simulated");
        // Sim_Hist_AB[i]->Draw("HIST SAME");

        Sim_Hist_conv_AB[i]->SetLineColor(kRed);
        Sim_Hist_conv_AB[i]->SetLineWidth(2);
        Sim_Hist_conv_AB[i]->SetTitle("Simulated Conv");
        Sim_Hist_conv_AB[i]->Draw("HIST SAME");    

        if (i == 0)
        {
            TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
            legend->AddEntry(Exp_Hist_calib_AB[i], "Experimental", "l");
            legend->AddEntry(Sim_Hist_conv_AB[i], "Simulated Conv", "l");
            legend->Draw("SAME");
        }    
    }
    
    c->Draw();
}