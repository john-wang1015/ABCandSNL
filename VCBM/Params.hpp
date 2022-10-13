#pragma once

#include <random>

std::uniform_real_distribution<double> uniform_random(0.0, 1.0);

class Params
{
private:
    std::random_device rd;
    std::mt19937 generator{ rd() };

public:
    //int created, opportunities;

    int gage, page;
    double p_0, dmax, p_psc;

    Params(double p_0, double p_psc, double dmax, int gage, int page)
    {
        this->p_0 = p_0;
        this->p_psc = p_psc;
        this->dmax = dmax;
        this->gage = gage;
        this->page = page;
    }

    double RandomDouble()
    {
        return uniform_random(generator);
    }

    bool WithProbability(double prob)
    {
        return RandomDouble() < prob;
    }

    // initialized below
    static const double s; // Cell-cell spring rest length
    static const double mu;
    static const double Delta_t;
    static const double rmin; // minimum distance needed between at least one cell and the cell proliferating
    static const int tinterval;
};

const double Params::s = 3.03;
const double Params::mu = 0.01;
const double Params::Delta_t = 30;
const double Params::rmin = 2;
const int Params::tinterval = 24;
