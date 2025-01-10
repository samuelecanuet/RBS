g++ -o RBS_All_forMCMC RBS_All_forMCMC.C `root-config --cflags --libs` 
g++ -o RBS_All_forMCMC_fit RBS_All_forMCMC_fit.C `root-config --cflags --libs` -lSpectrum 
