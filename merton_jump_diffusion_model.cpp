#include <iostream>
#include <random>
#include <vector>
using namespace std;

double gen_poisson(double lambda, mt19937 &gen)
{
    poisson_distribution<int> dist(lambda);
    double N = dist(gen);
    return N;
}

double gen_normal(int n, mt19937 &gen)
{
    normal_distribution<double> dist(0.0, 1.0);
    double z = dist(gen);
    return z;
}

vector<double> jump_diffusion_model(
    int n,
    double T,
    double short_rate,
    double volatility,
    double lambda,
    double jump_mean,
    double jump_volatility)
{
    random_device rd;
    mt19937 gen(rd());
    double dt = T / n;
    double k = exp(jump_mean + pow(jump_volatility, 2) / 2) - 1;
    vector<double> S(n + 1, 1.0);

    for (int i = 0; i < n; ++i)
    {
        double N = gen_poisson(lambda * dt, gen);
        double Z = gen_normal(n, gen);
        double X = gen_normal(N, gen);

        double jump = jump_mean * N + jump_volatility * sqrt(N) * X;
        S[i + 1] = S[i] * exp((short_rate - lambda * k - pow(volatility, 2) / 2) * dt + volatility * sqrt(dt) * Z + jump);
    }

    return S;
}