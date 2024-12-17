#include <iostream>
#include <vector>
#include <random>
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

// Gaussian function
double Gaussian(double x, double mean, double sigma) {
    return exp(-0.5 * pow((x - mean) / sigma, 2)) / (sigma * sqrt(2 * M_PI));
}

double f(double *x, double *par)

// Generate fake data
void GenerateFakeData(std::vector<double> &xData, std::vector<double> &yData, 
                      double mean, double sigma, int nPoints, double noiseLevel) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> gaussian(mean, sigma);
    std::uniform_real_distribution<> noise(-noiseLevel, noiseLevel);

    for (int i = 0; i < nPoints; ++i) {
        double x = -3.0 * sigma + i * (6.0 * sigma / nPoints);
        double y = Gaussian(x, mean, sigma) + noise(gen);
        xData.push_back(x);
        yData.push_back(y);
    }
}

int test() {
    // Generate fake data
    std::vector<double> xData, yData;
    double trueMean = 1.0;
    double trueSigma = 0.5;
    int nPoints = 50;
    double noiseLevel = 0.05;
    GenerateFakeData(xData, yData, trueMean, trueSigma, nPoints, noiseLevel);

    // Set up minimizer
    ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
    ROOT::Math::Functor functor(&Gaussian, 2);
    minimizer.SetMaxFunctionCalls(10000);
    minimizer.SetTolerance(1e-6);
    minimizer.SetPrintLevel(1);
    minimizer.SetFunction(functor);

    // Set initial parameter guesses
    minimizer.SetVariable(0, "mean", 1.0, 0.1);   // Initial guess for mean
    minimizer.SetVariable(1, "sigma", 0.5, 0.1);  // Initial guess for sigma

    // Perform minimization
    if (!minimizer.Minimize()) {
        std::cerr << "Minimization failed!" << std::endl;
        return 1;
    }

    // Get the results
    const double *params = minimizer.X();
    std::cout << "Fit results:" << std::endl;
    std::cout << "Mean = " << params[0] << std::endl;
    std::cout << "Sigma = " << params[1] << std::endl;

    // Plot the data and the fit
    TGraphErrors graph(nPoints);
    for (int i = 0; i < nPoints; ++i) {
        graph.SetPoint(i, xData[i], yData[i]);
    }
    graph.SetTitle("Gaussian Fit; x; y");
    graph.SetMarkerStyle(20);

    // Plot the Gaussian fit curve
    TCanvas canvas("canvas", "Gaussian Fit", 800, 600);
    graph.Draw("AP");

    TF1 fitFunction("fitFunction", "gaus", -3.0 * trueSigma, 3.0 * trueSigma);
    fitFunction.SetParameters(1.0, params[0], params[1]); // Set amplitude, mean, and sigma
    fitFunction.Draw("SAME");

    canvas.SaveAs("gaussian_fit.png");

    return 0;
}
