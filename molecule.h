#pragma once
#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>
#include <cassert>
using namespace std;
using namespace arma;

struct AO {
    string AO_type;
    string chemSym;
    int valence;
    rowvec center;
    ivec lmn;
    vec exponents;
    vec contraction_coeffs;
    vec norm_constants;

    AO(string AO_type, string chemSym, int valence, rowvec center, ivec lmn, vec exponents, vec contraction_coeffs);
    void calcNormConstants();
};

class Molecule {
    public:
        Molecule(const char* filename);
        int nAtoms;
        ivec atomicNumbers;
        vector<string> atomicSymbols;
        ivec atomValences;
        mat coordinates; 
        int charge;
        
        int hydrogenCount;
        int heavyAtomCount;
        int nElectrons; 
        int pAlpha;
        int qBeta;
        int nBasisFunctions;
        vector<AO> basisFunctionsList;
        map<int, int> aoIndexToAtom;

        map<string, vec> exponentsMap;
        map<string, vec> sContractionMap;
        map<string, vec> pContractionMap;

        void printMoleculeInfo();

    private:
        int countBasisFunctions();
        int countElectrons();
        vector<AO> buildBasisFunctionsList();
};

// Functions used to calculate overlap integrals
double overlapIntegral1D(double alpha, double beta, double center_a, double center_b, int lA, int lB);
double overlapIntegral3D(rowvec centers_a, rowvec centers_b, double alpha, double beta, ivec lmn_a, ivec lmn_b);
double calcContractedOverlap(AO AO1, AO AO2);
mat calcOverlapMatrix(Molecule molecule);
