#include <iostream>
#include <random>
#include <vector>
#include <cmath>
using namespace std;

double gen_gamma(double alpha, double beta, mt19937 &gen)
{
    gamma_distribution<double> dist(alpha, beta);
    return dist(gen);
}

double gen_normal(mt19937 &gen)
{
    normal_distribution<double> dist(0.0, 1.0);
    return dist(gen);
}

vector<double> variance_gamma_process(
    int n,
    double T,
    double short_rate,
    double sigma,
    double theta,
    double beta)
{
    random_device rd;
    mt19937 gen(rd());

    double dt = T / n;
    vector<double> S(n + 1, 1.0);
    double drift_correction = (1.0 / beta) * log(1.0 - theta * beta - 0.5 * pow(sigma, 2) * beta);

    for (int i = 0; i < n; ++i)
    {
        double G = gen_gamma(beta * dt, beta, gen);
        double Z = gen_normal(gen);

        double X = theta * G + sigma * sqrt(G) * Z;
        S[i + 1] = S[i] * exp((drift_correction + short_rate) * dt + X);
    }

    return S;
}
