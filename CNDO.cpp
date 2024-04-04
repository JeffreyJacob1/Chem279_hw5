#include "CNDO.h"
#include <cassert>
/**
 * CNDO/2 Hartree Fock method constructor
 * 
 */
CNDO::CNDO(Molecule molecule, mat overlapMatrix)
: molecule(molecule), overlapMatrix(overlapMatrix) {

    // Map the AO index to the atom it belongs to
    int index = 0;
    for (int A = 0; A < molecule.nAtoms; A++) {
        int numAOs = calcNumAOs(molecule.atomicSymbols[A]);
        for (int i = 0; i < numAOs; i++) {
            aoIndexToAtom[index] = A;
            index++;
        }
    }

    // Initialize semi-empirical parameters
    diagCNDOPara["H"]["1s"] = 7.176;

    diagCNDOPara["C"]["2s"] = 14.051;
    diagCNDOPara["C"]["2px"] = 5.572;
    diagCNDOPara["C"]["2py"] = 5.572;
    diagCNDOPara["C"]["2pz"] = 5.572;

    diagCNDOPara["N"]["2s"] = 19.316;
    diagCNDOPara["N"]["2px"] = 7.275;
    diagCNDOPara["N"]["2py"] = 7.275;
    diagCNDOPara["N"]["2pz"] = 7.275;

    diagCNDOPara["O"]["2s"] = 25.390;
    diagCNDOPara["O"]["2px"] = 9.111;
    diagCNDOPara["O"]["2py"] = 9.111;
    diagCNDOPara["O"]["2pz"] = 9.111;

    diagCNDOPara["F"]["2s"] = 32.272;
    diagCNDOPara["F"]["2px"] = 11.080;
    diagCNDOPara["F"]["2py"] = 11.080;
    diagCNDOPara["F"]["2pz"] = 11.080;

    offDiagCNDOPara["H"] = 9.;
    offDiagCNDOPara["C"] = 21.;
    offDiagCNDOPara["N"] = 25.;
    offDiagCNDOPara["O"] = 31.;
    offDiagCNDOPara["F"] = 39.;

    // Initialize matrices
    alphaCoeffMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    betaCoeffMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);

    alphaDensityMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    betaDensityMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    totalDensity = arma::zeros(molecule.nBasisFunctions);

    alphaEnergy.set_size(molecule.nBasisFunctions);
    betaEnergy.set_size(molecule.nBasisFunctions);

    gammaMatrix = calcGammaMatrix();
    alphaFockMat = calcFockMat(alphaDensityMat);
    betaFockMat = calcFockMat(betaDensityMat);
    hCoreMat = calcHCoreMat();

    nuclearRepulsionEnergy = calcNuclearRepulsionEnergy();
}

/**
 * 
 * 
 * This function calculates the number of AOs for an atom given the atomic
 * symbol. Hydrogen atoms have 1 AO, while all other atoms have 4 AOs.
 */
int CNDO::calcNumAOs(string chemSym) {
    if (chemSym == "H") {
        return 1;
    }
    else {
        return 4;
    }
}

/**
 * 
 * This function calculates the 2-electron repulsion integrals evaluated over the
 * square of the valence s orbital centered on atoms A and B.
 * 
 * param AO1 The first AO
 * param AO2 The second AO
 * 
 * return double The 2-electron repulsion integral
 */
double CNDO::calc2eIntegral(AO AO1, AO AO2) {
    if (!(accu(AO1.lmn) == 0 && accu(AO2.lmn) == 0)) {
        cout << "Error: 2e integrals only implemented for s orbitals" << endl;
    }

    // d prime = contraction coefficients * normalization constant
    vec dprime_a = AO1.contraction_coeffs % AO1.norm_constants;
    vec dprime_b = AO2.contraction_coeffs % AO2.norm_constants;

    int len = AO1.exponents.n_elem;

    double gamma = 0.0;
    for (int k1 = 0; k1 < len; k1++) {
        for (int k2 = 0; k2 < len; k2++) {
            double sigmaA = 1.0/(AO1.exponents(k1) + AO1.exponents(k2)); // eq 3.10
            for (int l1 = 0; l1 < len; l1++) {
                for (int l2 = 0; l2 < len; l2++) {
                    double sigmaB = 1.0/(AO2.exponents(l1) + AO2.exponents(l2)); // eq 3.10
                    
                    double I2e = pg2eIntegral(AO1.center, AO2.center, sigmaA, sigmaB);  // eq 3.14
                    gamma += dprime_a(k1) * dprime_a(k2) * dprime_b(l1) * dprime_b(l2) * I2e;
                }
            }
        }
    }
    return gamma;
}

/**
 * 
 * This function calculates the 2-electron integral over primitive Gaussians
 * given the centers, exponents, and contraction coefficients of the two Gaussians.
 * 
 * param center_a 
 * param center_b 
 * param sigmaA 
 * param sigmaB 
 * 
 * returns double 2-electron integral
 */
double CNDO::pg2eIntegral(rowvec center_a, rowvec center_b, double sigmaA, double sigmaB) {
    double U = pow(M_PI * sigmaA, 1.5) * pow(M_PI * sigmaB, 1.5);
    double V2 = 1.0 / (sigmaA + sigmaB);

    double distance = norm(center_a - center_b, 2);

    if (distance == 0.0) {
        return U * sqrt(2 * V2) * sqrt(2 / M_PI) * hartree2eV;  // eq 3.15
    } 

    double sqrtT = sqrt(V2) * distance;
    double result = U / distance * erf(sqrtT); 
    return result * hartree2eV;
}

/**
 * 
 * This function calculates the gamma matrix of 2-electron repulsion integrals
 * for the entire molecule using the calc2eIntegral function.
 * 
 * param molecule 
 * 
 * returns mat The gamma matrix
 */
mat CNDO::calcGammaMatrix() {
    // Make a list of only basis functions with s orbitals
    vector<AO> sBasisFunctionsList;
    for (int i = 0; i < molecule.nBasisFunctions; i++) {
        if (molecule.basisFunctionsList[i].AO_type == "1s" || molecule.basisFunctionsList[i].AO_type == "2s") {
            sBasisFunctionsList.push_back(molecule.basisFunctionsList[i]);
        }
    }

    mat gamma_matrix = arma::zeros<mat>(molecule.nAtoms, molecule.nAtoms);

    // Loop over all s orbital basis function combinations
    for (int i = 0; i < sBasisFunctionsList.size(); i++) {
        for (int j = 0; j < sBasisFunctionsList.size(); j++) {
            gamma_matrix(i, j) = calc2eIntegral(sBasisFunctionsList[i], sBasisFunctionsList[j]);
        }
    }
    return gamma_matrix;
}

/**
 * 
 * This function calculates the nuclear repulsion energy of the molecule
 * given the charges of the atoms
 * 
 * returns double The nuclear repulsion energy
 */
double CNDO::calcNuclearRepulsionEnergy() {
    double nuclearRepulsionEnergy = 0.0;
    for (int A = 0; A < molecule.nAtoms; A++) {
        for (int B = 0; B < A; B++) {
            double distance = norm(molecule.coordinates.row(A) - molecule.coordinates.row(B), 2);
            nuclearRepulsionEnergy += molecule.atomValences(A) * molecule.atomValences(B) / distance;
        }
    }
    return nuclearRepulsionEnergy * hartree2eV;
}

   
/**
 * 
 * This function calculates the total density vector using the alpha and beta
 * density matrices.
 * 
 * returns mat The total density vector
 */
vec CNDO::calcTotalDensity() {
    vec totalDensity = arma::zeros(molecule.nAtoms);

    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        // Get atom associated with mu AO
        int A = aoIndexToAtom[mu];
        totalDensity(A) += alphaDensityMat(mu, mu) + betaDensityMat(mu, mu);
    }

    return totalDensity;
}

/**
 * 
 * This function calculates the CNDO/2 Fock matrix given either the alpha or
 * beta density matrix.
 * 
 *  densityMat The density matrix (alpha or beta)
 * 
 * returns mat The CNDO/2 Fock matrix
 */
mat CNDO::calcFockMat(mat densityMat) {
    mat fockMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);

    // Loop over all AOs in molecule
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
            // Get atoms associated with mu and nu AOs
            int A = aoIndexToAtom[mu];
            int B = aoIndexToAtom[nu];
            string chemSymA = molecule.atomicSymbols[A];
            string chemSymB = molecule.atomicSymbols[B];
            double gammaAA = gammaMatrix(A, A);
            double gammaAB = gammaMatrix(A, B);
            double pAA = totalDensity(A);
            double ZA = molecule.atomValences(A);

            // Calculate the diagonal elements of the matrix
            if (mu == nu) {
                string AO_type = molecule.basisFunctionsList[mu].AO_type;
                fockMat(mu, nu) = -diagCNDOPara[chemSymA][AO_type] + \
                                  ((pAA - ZA) - (densityMat(mu, mu) - 0.5)) * gammaAA;
                
                // Update the diagonal elements of the matrix when A != B
                for (int B = 0; B < molecule.nAtoms; B++) {
                    if (A != B) {
                        double pBB = totalDensity(B);
                        double ZB = molecule.atomValences(B);
                        double gammaAB = gammaMatrix(A, B);
                        fockMat(mu, nu) += (pBB - ZB) * gammaAB;
                    }
                }
            }

            // Calculate the off-diagonal elements of the matrix
            else {
                fockMat(mu, nu) = (-offDiagCNDOPara[chemSymA] - offDiagCNDOPara[chemSymB]) \
                                  / 2.0 * overlapMatrix(mu, nu) - (densityMat(mu, nu) * gammaAB);
            }
        }
    }
    return fockMat;
}

/**
 * 
 * This function calculates the core Hamiltonian matrix similar to the Fock
 * matrix, but independent of electron density.
 * 
 * returns mat The core Hamiltonian matrix
 */
mat CNDO::calcHCoreMat() {
    mat hCoreMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    
    // Loop over all AOs in molecule
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
            int A = aoIndexToAtom[mu];
            int B = aoIndexToAtom[nu];
            string chemSymA = molecule.atomicSymbols[A];
            string chemSymB = molecule.atomicSymbols[B];
            double gammaAA = gammaMatrix(A, A);
            double gammaAB = gammaMatrix(A, B);
            double ZA = molecule.atomValences(A);

            // Calculate the diagonal elements of the matrix
            if (mu == nu) {
                string AO_type = molecule.basisFunctionsList[mu].AO_type;
                hCoreMat(mu, nu) = -diagCNDOPara[chemSymA][AO_type] - (ZA - 0.5) * gammaAA;

                for (int B = 0; B < molecule.nAtoms; B++) {
                    if (A != B) {
                        double ZB = molecule.atomValences(B);
                        double gammaAB = gammaMatrix(A, B);
                        hCoreMat(mu, nu) -= ZB * gammaAB;
                    }
                }
            }
            // Calculate the off-diagonal elements of the matrix
            else {
                hCoreMat(mu, nu) = (-offDiagCNDOPara[chemSymA] - offDiagCNDOPara[chemSymB]) \
                                   / 2.0 * overlapMatrix(mu, nu);
            }
        }
    }
    return hCoreMat;
}

/**
 * 
 * This function calculates the density matrix given the coefficient matrix
 * for either the alpha or beta electrons.
 * 
 *  coeffMat The coefficient matrix to use (alpha or beta)
 *  type The type of matrix (alpha or beta)
 * 
 * return mat The density matrix
 */
mat CNDO::calcDensityMat(mat coeffMatA, string type) {
    mat densityMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);

    if (type == "alpha") {
        for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
            for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
                for (int i = 0; i < molecule.pAlpha; i++) {
                    densityMat(mu, nu) += coeffMatA(mu, i) * coeffMatA(nu, i);
                }
            }
        }
    }

    else if (type == "beta") {
        for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
            for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
                for (int i = 0; i < molecule.qBeta; i++) {
                    densityMat(mu, nu) += coeffMatA(mu, i) * coeffMatA(nu, i);
                }
            }
        }
    }
    return densityMat;
}

/**
 * 
 * This function calculates the total energy after convergence is reached.
 * 
 * returns double The total energy
 */
double CNDO::calcTotalEnergy() {
    double totalEnergy = 0.0;
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
            totalEnergy += (alphaDensityMat(mu, nu) * (hCoreMat(mu, nu) + alphaFockMat(mu, nu)) + \
                            betaDensityMat(mu, nu) * (hCoreMat(mu, nu) + betaFockMat(mu, nu)));
        }
    }
    totalEnergy /= 2.0;
    totalEnergy += nuclearRepulsionEnergy;
    return totalEnergy;
}

/**
 * Self consistent field (SCF) cycle.
 * 
 * This function runs the SCF cycle until convergence is reached. The SCF cycle
 * 
 */
void CNDO::scfCycle() {
    // 1) Guess the density matrix
    alphaDensityMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    betaDensityMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);

    bool converged = false;
    int scfCycleCount = 0;  

    // 6) If not converged, repeat from step 2
    while (!converged) {
        scfCycleCount++;

        // 2) Calculate the Fock matrix
        alphaFockMat = calcFockMat(alphaDensityMat);
        betaFockMat = calcFockMat(betaDensityMat);

        // 3) Diagonalize the Fock matrix and obtain the coefficient matrix
        eig_sym(alphaEnergy, alphaCoeffMat, alphaFockMat);
        eig_sym(betaEnergy, betaCoeffMat, betaFockMat);

        // Save
        oldAlphaDensityMat = alphaDensityMat;
        oldBetaDensityMat = betaDensityMat;

        // 4) Calculate the density matrix
        alphaDensityMat = calcDensityMat(alphaCoeffMat, "alpha");
        betaDensityMat = calcDensityMat(betaCoeffMat, "beta");
        totalDensity = calcTotalDensity();

        // 5) Check for convergence (tolerance = 1e-6)
        if (abs(alphaDensityMat - oldAlphaDensityMat).max() < 1e-6 && \
            abs(betaDensityMat - oldBetaDensityMat).max() < 1e-6) {
    
            // 7) If converged, calculate the total energy
            converged = true;
            totalEnergy = calcTotalEnergy();
            cout << "Nuclear Repulsion Energy: " << nuclearRepulsionEnergy << " eV." << endl;
            cout << "Total Energy: " << totalEnergy << " eV." << endl;
        }
    }
}

/**
 * 
 * The expression for x involves isolating the terms that multiply the overlap matrix.
 * The derived expression for the X matrix is:
 * x(mu, nv) = (BetaA + BetaB) * density(mu, nv)
 * Where BetaA and BetaB refer to the off-diagonal parameters for atoms A and B.
 * 
 * return mat The X matrix
 */
mat CNDO::calcXMatrix() {
    mat X = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);

    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
            int A = aoIndexToAtom[mu];
            int B = aoIndexToAtom[nu];
            string chemSymA = molecule.atomicSymbols[A];
            string chemSymB = molecule.atomicSymbols[B];
            double betaA = offDiagCNDOPara[chemSymA];
            double betaB = offDiagCNDOPara[chemSymB];
            X(mu, nu) = (betaA + betaB) * (alphaDensityMat(mu, nu) + betaDensityMat(mu, nu));
        }
    }
    return X;
}

/**
 * 
 * The expression for y involves isolating the terms that multiply gammaAB.
 * The derived expression for the Y matrix is:
 * y(A, B) = totalDensity(A) * totalDensity(B) - ZB * totalDensity(A) - ZA * totalDensity(B)
 *           - sum(mu -> A) sum(nu -> B) (alphaDensity(mu, nu)**2 + betaDensity(mu, nu)**2)
 * 
 * returns mat The Y matrix
 */
mat CNDO::calcYMatrix() {
    mat Y = arma::zeros(molecule.nAtoms, molecule.nAtoms);

    for (int A = 0; A < molecule.nAtoms; A++) {
        for (int B = 0; B < molecule.nAtoms; B++) {
            double ZA = molecule.atomValences(A);
            double ZB = molecule.atomValences(B);

            Y(A, B) = totalDensity(A) * totalDensity(B) - ZB * totalDensity(A) - ZA * totalDensity(B); 
        }
    }

    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
            int A = aoIndexToAtom[mu];
            int B = aoIndexToAtom[nu];
            Y(A, B) -= (alphaDensityMat(mu, nu) * alphaDensityMat(mu, nu) + \
                        betaDensityMat(mu, nu) * betaDensityMat(mu, nu));
        }
    }

    return Y;
}

/**
 *
 * This function calculates the overlap integral gradient for one dimension and 
 * utilizes the overlapIntegral1D function from molecule.cpp to do so.
 * 
 * param center_a The center of the first primitive Gaussian (1D)
 * param center_b The center of the second primitive Gaussian (1D)
 * param alpha The exponent of the first primitive Gaussian 
 * param beta The exponent of the second primitive Gaussian
 * param lmn_a The angular momentum of the first primitive Gaussian
 * param lmn_b The angular momentum of the second primitive Gaussian
 * 
 * returns double The overlap integral gradient for one dimension
 */
double CNDO::calcOverlapGrad(double alpha, double beta, double center_a, double center_b, int lA, int lB) {
    if (lA == 0) {
        return 2 * alpha * overlapIntegral1D(alpha, beta, center_a, center_b, lA + 1, lB);
    }
    else {
        return -lA * overlapIntegral1D(alpha, beta, center_a, center_b, lA - 1, lB) + \
                2 * alpha * overlapIntegral1D(alpha, beta, center_a, center_b, lA + 1, lB);
    }
}

/**
 * 
 * This function calculates the overlap integral gradient for dx, dy, and dz for two primitive Gaussians. 
 * The expression for the overlap integral gradient is:
 * dS/dR = (dSx/dx * Sy * Sz, Sx * dSy/dy * Sz, Sx * Sy * dSz/dz)
 * 
 * param centers_a The center of the first primitive Gaussian
 * param centers_b The center of the second primitive Gaussian
 * param alpha The exponent of the first primitive Gaussian
 * param beta The exponent of the second primitive Gaussian
 * param lmn_a The angular momentum of the first primitive Gaussian
 * param lmn_b The angular momentum of the second primitive Gaussian
 * 
 * returns vec The overlap integral gradient for dx, dy, and dz
 */
vec CNDO::calcOverlapGradient3D(rowvec centers_a, rowvec centers_b, double alpha, double beta, ivec lmn_a, ivec lmn_b) {
    // Calculate the overlap integral for each dimension
    double Sx = overlapIntegral1D(alpha, beta, centers_a(0), centers_b(0), lmn_a(0), lmn_b(0));
    double Sy = overlapIntegral1D(alpha, beta, centers_a(1), centers_b(1), lmn_a(1), lmn_b(1));
    double Sz = overlapIntegral1D(alpha, beta, centers_a(2), centers_b(2), lmn_a(2), lmn_b(2));

    // Calculate the overlap integral gradient for each dimension
    double dSx = calcOverlapGrad(alpha, beta, centers_a(0), centers_b(0), lmn_a(0), lmn_b(0));
    double dSy = calcOverlapGrad(alpha, beta, centers_a(0), centers_b(0), lmn_a(1), lmn_b(1));
    double dSz = calcOverlapGrad(alpha, beta, centers_a(0), centers_b(0), lmn_a(2), lmn_b(2));

    vec overlap_grad = {dSx * Sy * Sz, Sx * dSy * Sz, Sx * Sy * dSz}; 
    return overlap_grad;
}

/**
 * 
 * This function includes the normalization constant and contraction coefficients in order to calculate
 * the derivative of the contracted overlap integral with respect to the nuclear center. The derivative 
 * is found for two given basis functions.
 * 
 *  AO1 The first AO
 *  AO2 The second AO
 * 
 * returns vec The derivative of the contracted overlap integral with respect to the nuclear center
 */
vec CNDO::calcContOverlapDeriv(AO AO1, AO AO2) {
    vec cont_overlap_grad = arma::zeros(3);

    // Look over all exponent, contraction coefficient, and normalization constant combinations
    for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
            vec unnorm_overlap_grad = calcOverlapGradient3D(AO1.center, AO2.center, AO1.exponents(k), 
                                                            AO2.exponents(l), AO1.lmn, AO2.lmn);
            cont_overlap_grad += AO1.contraction_coeffs(k) * AO2.contraction_coeffs(l) * \
                                 AO1.norm_constants(k) * AO2.norm_constants(l) * unnorm_overlap_grad;
        }
    }
    return cont_overlap_grad;
}

/**
 * 
 * This function calculates the overlap integral gradient matrix for the entire molecule.
 * The matrix is size (3, nBasisFunctions * nBasisFunctions) where the rows correspond to
 * the x, y, and z components of the gradient and the columns correspond to the basis
 * function combinations.
 * 
 * returns mat The overlap integral gradient matrix
 */
mat CNDO::calcOverlapGradientMat() {
    mat overlap_gradient_matrix = arma::zeros(3, molecule.nBasisFunctions * molecule.nBasisFunctions);

    // Loop over all basis function combinations
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
            AO AO1 = molecule.basisFunctionsList[mu];
            AO AO2 = molecule.basisFunctionsList[nu];
            
            vec cont_overlap_grad = calcContOverlapDeriv(AO1, AO2);
            overlap_gradient_matrix.col(mu * molecule.nBasisFunctions + nu) = cont_overlap_grad;
        }
    }
    return overlap_gradient_matrix;
}

/**
 * 
 * This function calculates the gamma 2e- integral gradient ([0]^RA) with respect to one dimension
 * of the nuclear center.
 * 
 *  center_a
 *  center_b
 *  sigmaA
 *  sigmaB
 * 
 * returns double 2e- integral gradient
 */
double CNDO::calcPG2eIntegralGrad(double center_a, double center_b, double sigmaA, double sigmaB) {
    double U = pow(M_PI * sigmaA, 1.5) * pow(M_PI * sigmaB, 1.5);
    double V = sqrt(1.0 / (sigmaA + sigmaB));

    double distance = abs(center_a - center_b);

    if (distance == 0.0) {
        return 0.0;
    }

    double T = V * V * distance * distance;
    double result = U / pow(distance, 2) * (2 * V / sqrt(M_PI) * exp(-T) - erf(sqrt(T))/distance) * (center_a - center_b);
    return result * hartree2eV;
}

/**
 *
 * This function calculates the gamma 2e- integral gradient for all dimensions (dx, dy, and dz) for two primitive Gaussians.
 * 
 * returns vec The gamma 2e- integral gradient for dx, dy, and dz
 */
vec CNDO::calcPG2eIntegralGrad3D(rowvec centers_a, rowvec centers_b, double sigmaA, double sigmaB) {
    // Calculate the gamma 2e integral for each dimension
    double gamma_x = calcPG2eIntegralGrad(centers_a(0), centers_b(0), sigmaA, sigmaB);
    double gamma_y = calcPG2eIntegralGrad(centers_a(1), centers_b(1), sigmaA, sigmaB);
    double gamma_z = calcPG2eIntegralGrad(centers_a(2), centers_b(2), sigmaA, sigmaB);

    vec gamma_grad = {gamma_x, gamma_y, gamma_z};
    return gamma_grad;
}

/**
 * 
 * This function calculates the gradient of the 2e- repulsion integral given two atomic orbital shells.
 * 
 *  AO1 The first AO
 *  AO2 The second AO
 * 
 * returns vec The gradient of the 2e- repulsion integral
 */
vec CNDO::calc2eIntegralGrad(AO AO1, AO AO2) {
    if (!(accu(AO1.lmn) == 0 && accu(AO2.lmn) == 0)) {
        cout << "Error: 2e integrals only implemented for s orbitals" << endl;
    }

    // d prime = contraction coefficients * normalization constant
    vec dprime_a = AO1.contraction_coeffs % AO1.norm_constants;
    vec dprime_b = AO2.contraction_coeffs % AO2.norm_constants;

    int len = AO1.exponents.n_elem;

    vec gamma_grad = arma::zeros(3);
    for (int k1 = 0; k1 < len; k1++) {
        for (int k2 = 0; k2 < len; k2++) {
            double sigmaA = 1.0/(AO1.exponents(k1) + AO1.exponents(k2)); // eq 3.10
            for (int l1 = 0; l1 < len; l1++) {
                for (int l2 = 0; l2 < len; l2++) {
                    double sigmaB = 1.0/(AO2.exponents(l1) + AO2.exponents(l2)); // eq 3.10
                    
                    vec I2e_grad_3D = calcPG2eIntegralGrad3D(AO1.center, AO2.center, sigmaA, sigmaB);  // eq 3.14
                    gamma_grad += dprime_a(k1) * dprime_a(k2) * dprime_b(l1) * dprime_b(l2) * I2e_grad_3D;
                }
            }
        }
    }
    return gamma_grad;
}

/**
 * This function calculates the gradient of the gamma matrix with respect to the nuclear center for the
 * entire molecule using the calc2eIntegralGrad function for each basis function combination. The gamma gradient
 * matrix size is (3, nAtoms * nAtoms) where the rows correspond to the x, y, and z components of the gradient
 * and the columns correspond to the atom combinations.
 * 
 */
mat CNDO::calcGammaGradientMat() {
    // Make a list of only basis functions with s orbitals
    vector<AO> sBasisFunctionsList;
    for (int i = 0; i < molecule.nBasisFunctions; i++) {
        if (molecule.basisFunctionsList[i].AO_type == "1s" || molecule.basisFunctionsList[i].AO_type == "2s") {
            sBasisFunctionsList.push_back(molecule.basisFunctionsList[i]);
        }
    }
    mat gamma_RA = arma::zeros(3, molecule.nAtoms * molecule.nAtoms);

    // Loop over all s orbital basis function combinations
    for (int i = 0; i < sBasisFunctionsList.size(); i++) {
        for (int j = 0; j < sBasisFunctionsList.size(); j++) {
            gamma_RA.col(i * molecule.nAtoms + j) = calc2eIntegralGrad(sBasisFunctionsList[i], sBasisFunctionsList[j]);
        }
    }
    return gamma_RA;
}

/**
 * 
 * This function calculates the nuclear energy gradient for two atoms given their
 * valence electrons and coordinates.
 * 
 * ZA The number of valence electrons for atom A
 *  ZB The number of valence electrons for atom B
 * RA The coordinates of atom A
 *  RB The coordinates of atom B
 */
vec CNDO::calcNuclearEnergyGrad(double ZA, double ZB, rowvec RA, rowvec RB) {
    vec nuclear_energy_grad = arma::zeros(3);
    
    for (int i = 0; i < 3; i++) {
        double diff = abs(RA(i) - RB(i));
        if (diff == 0.0) {
            nuclear_energy_grad(i) = 0.0;
        }
        else {
            nuclear_energy_grad(i) = (ZA * ZB) * -0.5 * pow(diff, -1.5) * 2 * diff;
        }
    }
    return nuclear_energy_grad * hartree2eV;
}

/**
 * 
 * This function calculates the nuclear energy gradient matrix for the entire molecule.
 * The matrix is size (3, nAtoms) where the rows correspond to the x, y, and z components
 * of the gradient and the columns correspond to each atom.
 * 
 * returns mat The nuclear energy gradient matrix
 */
mat CNDO::calcNuclearEnergyGradMat() {
    mat nuclear_energy_gradient_matrix = arma::zeros(3, molecule.nAtoms);

    // Loop over all atom combinations
    for (int A = 0; A < molecule.nAtoms; A++) {
        for (int B = 0; B < A; B++) {
            vec nuclear_energy_grad = calcNuclearEnergyGrad(molecule.atomValences(A), molecule.atomValences(B), \
                                                            molecule.coordinates.row(A), molecule.coordinates.row(B));
            nuclear_energy_gradient_matrix.col(A) += nuclear_energy_grad;
            nuclear_energy_gradient_matrix.col(B) -= nuclear_energy_grad;
        }
    }
    return nuclear_energy_gradient_matrix;
}

/**

 * 
 * This function calculates the gradient of the total energy with respect to the nuclear center.
 * The size of the gradient is (3, nAtoms) where the rows correspond to the x, y, and z components
 * 
 * returns: mat The gradient of the total energy
 */
mat CNDO::calcTotalEnergyGradMat() {
    mat energy_RA = arma::zeros(3, molecule.nAtoms);
    mat X = calcXMatrix();
    mat Y = calcYMatrix();
    mat overlap_RA = calcOverlapGradientMat();
    mat gamma_RA = calcGammaGradientMat();
    mat nuclear_energy_RA = calcNuclearEnergyGradMat();

    // Multiply X(mu,nu) and S_RA(mu,nu)
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
            for (int i = 0; i < 3; i++) {
                energy_RA(i, aoIndexToAtom[mu]) += X(mu, nu) * overlap_RA(i, mu * molecule.nBasisFunctions + nu);
            }
        }
    }

    // Multiply Y(A,B) and gamma_RA(A,B)
    for (int A = 0; A < molecule.nAtoms; A++) {
        for (int B = 0; B < molecule.nAtoms; B++) {
            if (A != B) {
                for (int i = 0; i < 3; i++) {
                    energy_RA(i, A) += Y(A, B) * gamma_RA(i, A * molecule.nAtoms + B);
                }
            }
        }
    }

    // Add nuclear energy gradient
    energy_RA += nuclear_energy_RA;

    return energy_RA;
}