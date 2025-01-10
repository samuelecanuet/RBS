#include <iostream>
#include <cmath>
#include <TMath.h>

// Function to compute the difference from the desired p using gammq
double chiSquareDiff(double delta, double p, int nu) {
    return TMath::Gamma(nu / 2.0, delta / 2.0) - p;
}

// Bisection method for root-finding
double findChiSquareQuantile(double p, int nu, double tol = 1e-8, int max_iter = 1e8) {
    double a = 0;        // Lower bound
    double b = 20;     // Upper bound, adjust if necessary
    double c = a + (b - a) / 2.0;
    int iter = 0;

    while (iter < max_iter) {
        double fa = chiSquareDiff(a, p, nu);
        double fc = chiSquareDiff(c, p, nu);

        if (std::fabs(fc) < tol) {
            return c;  // c is the root
        }

        // Bisection method: check which subinterval to continue
        if (fa * fc < 0) {
            b = c;
        } else {
            a = c;
        }

        c = a + (b - a) / 2.0;
        iter++;
    }

    std::cerr << "Max iterations reached, result may not be accurate!" << std::endl;
    return c;  // Return the last estimate
}

int testgammaq() {
    double p = 0.682689492137086;  // Desired probability
    int nu = 8;      // Degrees of freedom

    double delta = findChiSquareQuantile(p, nu);
    std::cout << "For p = " << p << " and nu = " << nu << ", the chi-square quantile is: " << delta << std::endl;

    return 0;
}
