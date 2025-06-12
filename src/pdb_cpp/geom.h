#ifndef GEOM_H
#define GEOM_H

#pragma once
#include <cstring>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>

class CrystalPack {

public:

    void set_CRYST1_pdb(const std::string& line) {
        // Parse CRYST1 line
        alpha = std::stof(line.substr(6, 9));
        beta = std::stof(line.substr(15, 9));
        gamma = std::stof(line.substr(24, 9));
        a = std::stof(line.substr(33, 7));
        b = std::stof(line.substr(40, 7));
        c = std::stof(line.substr(47, 7));
        sGroup = line.substr(56, 10);
        try
        {
            z = std::stoi(line.substr(67, 3));
        }
        catch(const std::exception& e)
        {
            z = 1;
        }
        
    }

    std::string get_pdb_crystal_pack() const {
        std::stringstream ss;
        ss << "CRYST1"
           << std::setw(9) << std::setprecision(3) << std::fixed << alpha
           << std::setw(9) << std::setprecision(3) << std::fixed << beta
           << std::setw(9) << std::setprecision(3) << std::fixed << gamma
           << std::setw(7) << std::setprecision(2) << std::fixed << a
           << std::setw(7) << std::setprecision(2) << std::fixed << b
           << std::setw(7) << std::setprecision(2) << std::fixed << c
           << " P" << sGroup << "  "
           << std::setw(2) << z << "\n";
        return ss.str();
    }


    void clear() {
        alpha = beta = gamma = a = b = c = nan("");
    }

private:

    float alpha=nan(""), beta=nan(""), gamma=nan(""), a=nan(""), b=nan(""), c=nan("");
    int z;
    std::string sGroup;

};


class Transformation {
public:
    void parse_pdb_transformation(const std::string& text) {
        chains.clear();
        matrix.clear();
        std::istringstream iss(text);
        std::string line;

        while (std::getline(iss, line)) {
            if (line.size() >= 42 && line.substr(34, 7) == "CHAINS:") {
                std::string chains_str = line.substr(42);
                std::istringstream chainss(chains_str);
                std::string chain;
                while (std::getline(chainss, chain, ',')) {
                    // Remove leading/trailing whitespace
                    chain.erase(chain.begin(), std::find_if(chain.begin(), chain.end(), [](unsigned char ch) { return !std::isspace(ch); }));
                    chain.erase(std::find_if(chain.rbegin(), chain.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), chain.end());
                    if (!chain.empty())
                        chains.push_back(chain);
                }
            } else if (line.size() >= 19 && line.substr(0, 18) == "REMARK 350   BIOMT") {
                std::istringstream vals(line.substr(19));
                std::vector<float> row;
                float val;
                while (vals >> val) {
                    row.push_back(val);
                }
                matrix.push_back(row);
            }
        }
    }
    void print() const {
        std::cout << "Chains: ";
        for (const auto& chain : chains) {
            std::cout << chain << " ";
        }
        std::cout << std::endl;

        std::cout << "Matrix:" << std::endl;
        for (const auto& row : matrix) {
            for (const auto& val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }

    void clear() {
        chains.clear();
        matrix.clear();
    }
private:
    std::vector<std::string> chains;
    std::vector<std::vector<float>> matrix;
};

class Symmetry {
public:
    void parse_pdb_symmetry(const std::string& text) {
        matrix.clear();
        std::istringstream iss(text);
        std::string line;
        while (std::getline(iss, line)) {
            if (line.size() >= 19 && line.substr(0, 18) == "REMARK 290   SMTRY") {
                std::istringstream vals(line.substr(19));
                std::vector<float> row;
                float val;
                while (vals >> val) {
                    row.push_back(val);
                }
                matrix.push_back(row);
            }
        }
    }
    void print() const {
        std::cout << "Symmetry Matrix:" << std::endl;
        for (const auto& row : matrix) {
            for (const auto& val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }
    void clear() {
        matrix.clear();
    }
private:
    std::vector<std::vector<float>> matrix;
};

inline float calculate_distance(float x1, float y1, float z1, float x2, float y2, float z2){
    return std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}
inline float calculate_square_distance(float x1, float y1, float z1, float x2, float y2, float z2) {
    return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
}

/**************************************************************************
Implementation of Kabsch algorithm for finding the best rotation matrix
---------------------------------------------------------------------------
x    - x(i,m) are coordinates of atom m in set x            (input)
y    - y(i,m) are coordinates of atom m in set y            (input)
n    - n is number of atom pairs                            (input)
mode  - 0:calculate rms only                                (input)
        1:calculate u,t only                                (takes medium)
        2:calculate rms,u,t                                 (takes longer)
rms   - sum of w*(ux+t-y)**2 over all atom pairs            (output)
u    - u(i,j) is   rotation  matrix for best superposition  (output)
t    - t(i)   is translation vector for best superposition  (output)
**************************************************************************/
inline bool Kabsch(double **x, double **y, int n, int mode, double *rms,
    double t[3], double u[3][3]) {
    int i, j, m, m1, l, k;
    double e0, rms1, d, h, g;
    double cth, sth, sqrth, p, det, sigma;
    double xc[3], yc[3];
    double a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
    double sqrt3 = 1.73205080756888, tol = 0.01;
    int ip[] = { 0, 1, 3, 1, 2, 4, 3, 4, 5 };
    int ip2312[] = { 1, 2, 0, 1 };

    int a_failed = 0, b_failed = 0;
    double epsilon = 0.00000001;

    //initialization
    *rms = 0;
    rms1 = 0;
    e0 = 0;
    double c1[3], c2[3];
    double s1[3], s2[3];
    double sx[3], sy[3], sz[3];
    for (i = 0; i < 3; i++) {
        s1[i] = 0.0;
        s2[i] = 0.0;

        sx[i] = 0.0;
        sy[i] = 0.0;
        sz[i] = 0.0;
    }

    for (i = 0; i<3; i++) {
        xc[i] = 0.0;
        yc[i] = 0.0;
        t[i] = 0.0;
        for (j = 0; j<3; j++) {
            u[i][j] = 0.0;
            r[i][j] = 0.0;
            a[i][j] = 0.0;
            if (i == j) {
                u[i][j] = 1.0;
                a[i][j] = 1.0;
            }
        }
    }

    if (n<1) return false;

    //compute centers for vector sets x, y
    for (i = 0; i<n; i++) {
        for (j = 0; j < 3; j++) {
            c1[j] = x[i][j];
            c2[j] = y[i][j];

            s1[j] += c1[j];
            s2[j] += c2[j];
        }

        for (j = 0; j < 3; j++) {
            sx[j] += c1[0] * c2[j];
            sy[j] += c1[1] * c2[j];
            sz[j] += c1[2] * c2[j];
        }
    }
    for (i = 0; i < 3; i++) {
        xc[i] = s1[i] / n;
        yc[i] = s2[i] / n;
    }
    if (mode == 2 || mode == 0)
        for (int mm = 0; mm < n; mm++)
            for (int nn = 0; nn < 3; nn++)
                e0 += (x[mm][nn] - xc[nn]) * (x[mm][nn] - xc[nn]) + 
                      (y[mm][nn] - yc[nn]) * (y[mm][nn] - yc[nn]);
    for (j = 0; j < 3; j++){
        r[j][0] = sx[j] - s1[0] * s2[j] / n;
        r[j][1] = sy[j] - s1[1] * s2[j] / n;
        r[j][2] = sz[j] - s1[2] * s2[j] / n;
    }

    //compute determinant of matrix r
    det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])\
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])\
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    sigma = det;

    //compute tras(r)*r
    m = 0;
    for (j = 0; j<3; j++){
        for (i = 0; i <= j; i++){
            rr[m] = r[0][i] * r[0][j] + r[1][i] * r[1][j] + r[2][i] * r[2][j];
            m++;
        }
    }

    double spur = (rr[0] + rr[2] + rr[5]) / 3.0;
    double cof = (((((rr[2] * rr[5] - rr[4] * rr[4]) + rr[0] * rr[5])\
        - rr[3] * rr[3]) + rr[0] * rr[2]) - rr[1] * rr[1]) / 3.0;
    det = det*det;

    for (i = 0; i<3; i++) e[i] = spur;

    if (spur>0){
        d = spur*spur;
        h = d - cof;
        g = (spur*cof - det) / 2.0 - spur*h;

        if (h>0) {
            sqrth = sqrt(h);
            d = h*h*h - g*g;
            if (d<0.0) d = 0.0;
            d = atan2(sqrt(d), -g) / 3.0;
            cth = sqrth * cos(d);
            sth = sqrth*sqrt3*sin(d);
            e[0] = (spur + cth) + cth;
            e[1] = (spur - cth) + sth;
            e[2] = (spur - cth) - sth;

            if (mode != 0) {//compute a                
                for (l = 0; l<3; l = l + 2) {
                    d = e[l];
                    ss[0] = (d - rr[2]) * (d - rr[5]) - rr[4] * rr[4];
                    ss[1] = (d - rr[5]) * rr[1] + rr[3] * rr[4];
                    ss[2] = (d - rr[0]) * (d - rr[5]) - rr[3] * rr[3];
                    ss[3] = (d - rr[2]) * rr[3] + rr[1] * rr[4];
                    ss[4] = (d - rr[0]) * rr[4] + rr[1] * rr[3];
                    ss[5] = (d - rr[0]) * (d - rr[2]) - rr[1] * rr[1];

                    if (fabs(ss[0]) <= epsilon) ss[0] = 0.0;
                    if (fabs(ss[1]) <= epsilon) ss[1] = 0.0;
                    if (fabs(ss[2]) <= epsilon) ss[2] = 0.0;
                    if (fabs(ss[3]) <= epsilon) ss[3] = 0.0;
                    if (fabs(ss[4]) <= epsilon) ss[4] = 0.0;
                    if (fabs(ss[5]) <= epsilon) ss[5] = 0.0;

                    if (fabs(ss[0]) >= fabs(ss[2])) {
                        j = 0;
                        if (fabs(ss[0]) < fabs(ss[5])) j = 2;
                    }
                    else if (fabs(ss[2]) >= fabs(ss[5])) j = 1;
                    else j = 2;

                    d = 0.0;
                    j = 3 * j;
                    for (i = 0; i<3; i++) {
                        k = ip[i + j];
                        a[i][l] = ss[k];
                        d = d + ss[k] * ss[k];
                    }


                    //if( d > 0.0 ) d = 1.0 / sqrt(d);
                    if (d > epsilon) d = 1.0 / sqrt(d);
                    else d = 0.0;
                    for (i = 0; i<3; i++) a[i][l] = a[i][l] * d;
                }//for l

                d = a[0][0] * a[0][2] + a[1][0] * a[1][2] + a[2][0] * a[2][2];
                if ((e[0] - e[1]) >(e[1] - e[2])) {
                    m1 = 2;
                    m = 0;
                } else {
                    m1 = 0;
                    m = 2;
                }
                p = 0;
                for (i = 0; i<3; i++) {
                    a[i][m1] = a[i][m1] - d*a[i][m];
                    p = p + a[i][m1] * a[i][m1];
                }
                if (p <= tol) {
                    p = 1.0;
                    for (i = 0; i<3; i++) {
                        if (p < fabs(a[i][m])) continue;
                        p = fabs(a[i][m]);
                        j = i;
                    }
                    k = ip2312[j];
                    l = ip2312[j + 1];
                    p = sqrt(a[k][m] * a[k][m] + a[l][m] * a[l][m]);
                    if (p > tol) {
                        a[j][m1] = 0.0;
                        a[k][m1] = -a[l][m] / p;
                        a[l][m1] = a[k][m] / p;
                    }
                    else a_failed = 1;
                }//if p<=tol
                else {
                    p = 1.0 / sqrt(p);
                    for (i = 0; i<3; i++) a[i][m1] = a[i][m1] * p;
                }//else p<=tol  
                if (a_failed != 1) {
                    a[0][1] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
                    a[1][1] = a[2][2] * a[0][0] - a[2][0] * a[0][2];
                    a[2][1] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
                }
            }//if(mode!=0)       
        }//h>0

        //compute b anyway
        if (mode != 0 && a_failed != 1)//a is computed correctly
        {
            //compute b
            for (l = 0; l<2; l++) {
                d = 0.0;
                for (i = 0; i<3; i++) {
                    b[i][l] = r[i][0] * a[0][l] + 
                              r[i][1] * a[1][l] + r[i][2] * a[2][l];
                    d = d + b[i][l] * b[i][l];
                }
                //if( d > 0 ) d = 1.0 / sqrt(d);
                if (d > epsilon) d = 1.0 / sqrt(d);
                else d = 0.0;
                for (i = 0; i<3; i++) b[i][l] = b[i][l] * d;
            }
            d = b[0][0] * b[0][1] + b[1][0] * b[1][1] + b[2][0] * b[2][1];
            p = 0.0;

            for (i = 0; i<3; i++) {
                b[i][1] = b[i][1] - d*b[i][0];
                p += b[i][1] * b[i][1];
            }

            if (p <= tol) {
                p = 1.0;
                for (i = 0; i<3; i++) {
                    if (p<fabs(b[i][0])) continue;
                    p = fabs(b[i][0]);
                    j = i;
                }
                k = ip2312[j];
                l = ip2312[j + 1];
                p = sqrt(b[k][0] * b[k][0] + b[l][0] * b[l][0]);
                if (p > tol) {
                    b[j][1] = 0.0;
                    b[k][1] = -b[l][0] / p;
                    b[l][1] = b[k][0] / p;
                }
                else b_failed = 1;
            }//if( p <= tol )
            else {
                p = 1.0 / sqrt(p);
                for (i = 0; i<3; i++) b[i][1] = b[i][1] * p;
            }
            if (b_failed != 1) {
                b[0][2] = b[1][0] * b[2][1] - b[1][1] * b[2][0];
                b[1][2] = b[2][0] * b[0][1] - b[2][1] * b[0][0];
                b[2][2] = b[0][0] * b[1][1] - b[0][1] * b[1][0];
                //compute u
                for (i = 0; i<3; i++)
                    for (j = 0; j<3; j++)
                        u[i][j] = b[i][0] * a[j][0] + 
                                  b[i][1] * a[j][1] + b[i][2] * a[j][2];
            }

            //compute t
            for (i = 0; i<3; i++)
                t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1]) - 
                                                    u[i][2] * xc[2];
        }//if(mode!=0 && a_failed!=1)
    }//spur>0
    else //just compute t and errors
    {
        //compute t
        for (i = 0; i<3; i++)
            t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1]) - 
                                                u[i][2] * xc[2];
    }//else spur>0 

    //compute rms
    for (i = 0; i<3; i++) {
        if (e[i] < 0) e[i] = 0;
        e[i] = sqrt(e[i]);
    }
    d = e[2];
    if (sigma < 0.0) d = -d;
    d = (d + e[1]) + e[0];

    if (mode == 2 || mode == 0) {
        rms1 = (e0 - d) - d;
        if (rms1 < 0.0) rms1 = 0.0;
    }

    *rms = rms1;
    return true;
}

// Wrapper function to use Kabsch algorithm with std::vector and std::array
// This maintains compatibility with the existing quaternion_rotate interface
inline std::array<std::array<float, 3>, 3> quaternion_rotate(
    const std::vector<std::array<float, 3>>& X, 
    const std::vector<std::array<float, 3>>& Y) {
    
    if (X.size() != Y.size()) {
        throw std::runtime_error("Input arrays X and Y must have the same length");
    }
    
    size_t n = X.size();
    if (n == 0) {
        // Return identity matrix for empty input
        return {{{{1.0f, 0.0f, 0.0f}}, {{0.0f, 1.0f, 0.0f}}, {{0.0f, 0.0f, 1.0f}}}};
    }
    
    // Allocate memory for Kabsch algorithm
    double **x = new double*[n];
    double **y = new double*[n];
    for (size_t i = 0; i < n; i++) {
        x[i] = new double[3];
        y[i] = new double[3];
        for (int j = 0; j < 3; j++) {
            x[i][j] = static_cast<double>(X[i][j]);
            y[i][j] = static_cast<double>(Y[i][j]);
        }
    }
    
    double rms;
    double t[3];
    double u[3][3];
    
    // Call Kabsch algorithm (mode=1 to calculate rotation matrix only)
    bool success = Kabsch(x, y, static_cast<int>(n), 1, &rms, t, u);
    
    // Clean up memory
    for (size_t i = 0; i < n; i++) {
        delete[] x[i];
        delete[] y[i];
    }
    delete[] x;
    delete[] y;
    
    if (!success) {
        // Return identity matrix if Kabsch fails
        return {{{{1.0f, 0.0f, 0.0f}}, {{0.0f, 1.0f, 0.0f}}, {{0.0f, 0.0f, 1.0f}}}};
    }
    
    // Convert double matrix to float array
    std::array<std::array<float, 3>, 3> rotation_matrix;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rotation_matrix[i][j] = static_cast<float>(u[i][j]);
        }
    }
    
    return rotation_matrix;
}

#endif // GEOM_H