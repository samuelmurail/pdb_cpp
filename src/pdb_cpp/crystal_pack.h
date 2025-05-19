#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <array>

using namespace std;

class CrystalPack {

public:

    void setAlpha(float alpha) { this->alpha = alpha; }
    void setBeta(float beta) { this->beta = beta; }
    void setGamma(float gamma) { this->gamma = gamma; }
    void setA(float a) { this->a = a; }
    void setB(float b) { this->b = b; }
    void setC(float c) { this->c = c; }

    void clear() {
        alpha = 0.0;
        beta = 0.0;
        gamma = 0.0;
        a = 0.0;
        b = 0.0;
        c = 0.0;
    }

private:

    float alpha, beta, gamma, a, b, c;

    
};
