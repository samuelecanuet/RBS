void root2txt()
{
    TFile *f = new TFile("RBS_Results.root", "READ");

    TH1D* h = (TH1D*)f->Get("Writing/Sim_Hist_conv_1.2_85.0_550.0_A_-2_4keV");

    ofstream file;
    string filename = "SIM.txt";
    file.open(filename.c_str());  
    for (int bin = 0; bin < h->GetNbinsX(); bin++)
    {
        file << h->GetBinContent(bin) << endl;
    }

    file.close();
}