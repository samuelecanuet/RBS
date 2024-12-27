#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include <vector>
#include <cmath>

// Define the model function
double modelFunc(double *x, double *par) {
    return par[0] + par[1] * x[0] + par[2] * std::sin(par[3] * x[0]);
}

// Compute the chi-squared value
double computeChi2(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& errors, const double* params) {
    double chi2 = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        double *xPtr = const_cast<double*>(&x[i]);
        double diff = y[i] - modelFunc(xPtr, const_cast<double*>(params));
        chi2 += (diff * diff) / (errors[i] * errors[i]);
    }
    return chi2;
}

// MCMC sampling using Metropolis-Hastings
void test() {
    // Data points
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {1.1, 2.0, 2.9, 4.1, 5.2};
    std::vector<double> errors = {0.1, 0.1, 0.1, 0.1, 0.1};

    // Initial fit parameters and ranges
    double params[4] = {1.0, 1.0, 0.0, 1.0}; // Initial guess
    double stepSize[4] = {0.05, 0.05, 0.01, 0.05}; // Step size for MCMC
    double accepted = 0, total = 0;

    // Random number generator
    TRandom3 rand(0);

    // Store MCMC samples
    const int nSteps = 10000;
    std::vector<std::vector<double>> samples(nSteps, std::vector<double>(4));

    // Run MCMC
    for (int step = 0; step < nSteps; ++step) {
        // Propose new parameters
        double newParams[4];
        for (int i = 0; i < 4; ++i) {
            newParams[i] = params[i] + stepSize[i] * rand.Gaus();
        }

        // Calculate chi-squared for old and new parameters
        double chi2Old = computeChi2(x, y, errors, params);
        double chi2New = computeChi2(x, y, errors, newParams);

        // Accept or reject
        double alpha = std::exp(-0.5 * (chi2New - chi2Old));
        if (alpha >= 1.0 || rand.Uniform() < alpha) {
            for (int i = 0; i < 4; ++i) {
                params[i] = newParams[i];
            }
            accepted++;
        }
        total++;

        // Store sample
        for (int i = 0; i < 4; ++i) {
            samples[step][i] = params[i];
        }
    }

    // Analyze MCMC samples
    std::vector<double> means(4, 0.0);
    std::vector<double> variances(4, 0.0);
    for (const auto& sample : samples) {
        for (int i = 0; i < 4; ++i) {
            means[i] += sample[i];
            variances[i] += sample[i] * sample[i];
        }
    }
    for (int i = 0; i < 4; ++i) {
        means[i] /= nSteps;
        variances[i] = std::sqrt(variances[i] / nSteps - means[i] * means[i]);
    }

    // Output results
    for (int i = 0; i < 4; ++i) {
        std::cout << "Parameter " << i << ": Mean = " << means[i] << ", StdDev = " << variances[i] << std::endl;
    }
    std::cout << "Acceptance Ratio: " << accepted / total << std::endl;
}
