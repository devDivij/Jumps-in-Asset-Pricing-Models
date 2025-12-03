#include <iostream>
#include <random>
#include <vector>
#include <cmath>

using namespace std;

double gen_normal(mt19937 &gen)
{
    normal_distribution<double> dist(0.0, 1.0);
    return dist(gen);
}

double gen_uniform(mt19937 &gen)
{
    uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(gen);
}

double gen_inverse_gaussian(double delta, double gamma, mt19937 &gen)
{

    double Z = gen_normal(gen);
    double V = Z * Z;
    double a = 1.0 / gamma;

    double b = a * delta;
    double xi = a * V;
    double Y = a * (delta + (xi / 2.0) + sqrt(xi * (delta + (xi / 4.0))));

    double p = delta / (delta + gamma * Y);
    double U = gen_uniform(gen);

    if (U > p)
    {
        Y = (b * b) / Y;
    }

    return Y;
}

vector<double> normal_inverse_gaussian_process(
    int n,
    double T,
    double short_rate,
    double delta,
    double gamma,
    double beta)
{
    random_device rd;
    mt19937 gen(rd());

    double dt = T / n;
    vector<double> S(n + 1, 1.0);
    double drift_correction = -delta * (gamma - sqrt(gamma * gamma - 2.0 * beta - 1.0));

    for (int i = 0; i < n; ++i)
    {
        double Y = gen_inverse_gaussian(delta * dt, gamma, gen);
        double Z = gen_normal(gen);
        double X = beta * Y + sqrt(Y) * Z;
        S[i + 1] = S[i] * exp((short_rate + drift_correction) * dt + X);
    }

    return S;
}
