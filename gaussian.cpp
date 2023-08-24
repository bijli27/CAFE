#include "gaussian.h"

double guassianclass::calculateCumulativeProbability(double x, double mu, double sigma) {
 double z = (x - mu) / (sigma * std::sqrt(2)); // Calculate the z-score

    double cumulativeProbability = 0.5 * (1 + std::erf(z / std::sqrt(2))); // Calculate the cumulative probability
    
    return cumulativeProbability;
}


