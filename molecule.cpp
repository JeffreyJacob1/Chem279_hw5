#include "molecule.h"
#include "factorial.h"

/**
 * AO struct constructor
 * 
 */
AO::AO(string AO_type, string chemSym, int valence, rowvec center, ivec lmn, vec exponents, vec contraction_coeffs)
    : AO_type(AO_type), chemSym(chemSym), valence(valence), center(center), lmn(lmn), exponents(exponents), 
      contraction_coeffs(contraction_coeffs) {
    // Calculate the normalization constants
    calcNormConstants();
}

/**
 * Calculate the three normalization constants for a given atomic orbital (AO).
 * 
 */
void AO::calcNormConstants() {
    // Calculate the overlap integrals for each exponent with itself
    double selfOverlap_1 = overlapIntegral3D(center, center, exponents(0), exponents(0), lmn, lmn);
    double selfOverlap_2 = overlapIntegral3D(center, center, exponents(1), exponents(1), lmn, lmn);
    double selfOverlap_3 = overlapIntegral3D(center, center, exponents(2), exponents(2), lmn, lmn);  

    // Calculate the normalization constants
    norm_constants = {1.0 / sqrt(selfOverlap_1), 1.0 / sqrt(selfOverlap_2), 1.0 / sqrt(selfOverlap_3)};
}


/**
 *Molecule class constructor
 * 
 * 
 */
Molecule::Molecule(const char* filename) {
    // Open the file
    ifstream infile(filename);
    assert(infile.good());

    // Read in the number of atoms and charge
    infile >> nAtoms >> charge;

    atomicNumbers.resize(nAtoms);
    atomicSymbols.resize(nAtoms);
    atomValences.resize(nAtoms);
    coordinates.resize(nAtoms, 3);

    string line; 
    getline(infile, line);

    for (int i = 0; i < nAtoms; i++) {
        getline(infile, line);
        istringstream iss(line);

        // Read and store the atomic number
        int atomicNumber;
        iss >> atomicNumber;
        atomicNumbers(i) = atomicNumber;

        // Loop over the three coordinates x, y, z
        for (int j = 0; j < 3; j++) {
            iss >> coordinates(i, j);
        }
    }

    // Close the file
    infile.close();

    // Set the atomic symbols and valence electrons given the atomic numbers
    map<int, string> atomicSymbolsMap = {{1, "H"},
                                         {6, "C"},
                                         {7, "N"},
                                         {8, "O"},
                                         {9, "F"}};

    map<int, int> valenceElectronsMap = {{1, 1},
                                         {6, 4},
                                         {7, 5},
                                         {8, 6},
                                         {9, 7}};

    for (int i = 0; i < nAtoms; i++) {
        atomicSymbols[i] = atomicSymbolsMap[atomicNumbers(i)];
        atomValences(i) = valenceElectronsMap[atomicNumbers(i)];
    }

    // Count the heavy atoms and hydrogens
    hydrogenCount = 0;
    heavyAtomCount = 0;
    for (int i = 0; i < nAtoms; i++) {
        if (atomicNumbers(i) == 1) {
            hydrogenCount++;
        } else {
            heavyAtomCount++;
        }
    }

    // Calculate the number of electrons and basis functions
    nBasisFunctions = countBasisFunctions();
    nElectrons = countElectrons();
    int remainder = nElectrons % 2;
    pAlpha = nElectrons / 2 + remainder;
    qBeta = nElectrons / 2;

    // Exponents and contraction coefficients given
    exponentsMap["H"] = {3.42525091, 0.62391373, 0.16885540};
    sContractionMap["H"] = {0.15432897, 0.53532814, 0.44463454};

    exponentsMap["C"] = {2.94124940, 0.68348310, 0.22228990};
    sContractionMap["C"] = {-0.09996723, 0.39951283, 0.70011547};
    pContractionMap["C"] = {0.15591627, 0.60768372, 0.39195739};

    exponentsMap["N"] = {3.78045590, 0.87849660, 0.28571440};
    sContractionMap["N"] = {-0.09996723, 0.39951283, 0.70011547};
    pContractionMap["N"] = {0.15591627, 0.60768372, 0.39195739};

    exponentsMap["O"] = {5.03315130, 1.16959610, 0.38038900};
    sContractionMap["O"] = {-0.09996723, 0.39951283, 0.70011547};
    pContractionMap["O"] = {0.15591627, 0.60768372, 0.39195739};

    exponentsMap["F"] = {6.46480320, 1.50228120, 0.48858850};
    sContractionMap["F"] = {-0.09996723, 0.39951283, 0.70011547};
    pContractionMap["F"] = {0.15591627, 0.60768372, 0.39195739};

    // Build the list of basic functions
    basisFunctionsList = buildBasisFunctionsList();
}

/**
 *  Evaluate the number of basis functions in the molecule
 */
int Molecule::countBasisFunctions() {
    return 4 * heavyAtomCount + hydrogenCount;
}

/**
 *  Evaluate the number of valence electrons in the moleculex
 */
int Molecule::countElectrons() {
    int nElectrons = 0;

    for (int i = 0; i < nAtoms; i++) {
        nElectrons += atomValences(i);
    }
    return nElectrons - charge;
}

/**
 *  Create a list of all the basis functions in the molecule, which are contracted Gaussians.
 *
 * return vector<AO> A vector of basic functions
 */
vector<AO> Molecule::buildBasisFunctionsList() {
    vector<AO> basisFunctionsList;

    for (int i = 0; i < nAtoms; i++) {
        string atomSym = atomicSymbols[i];
        int valence = atomValences(i);
        string orbType;
        if (atomSym == "H") {
            orbType = "1s";
        } else {
            orbType = "2s";
        }

        AO AO_s = {orbType, atomSym, valence, coordinates.row(i), {0, 0, 0}, exponentsMap[atomSym], sContractionMap[atomSym]};
        basisFunctionsList.push_back(AO_s);

        if (orbType == "2s") {
            // Add p-type basis functions for elements with "2s" contraction
            AO AO_2px = {"2px", atomSym, valence, coordinates.row(i), {1, 0, 0}, exponentsMap[atomSym], pContractionMap[atomSym]};
            AO AO_2py = {"2py", atomSym, valence, coordinates.row(i), {0, 1, 0}, exponentsMap[atomSym], pContractionMap[atomSym]};
            AO AO_2pz = {"2pz", atomSym, valence, coordinates.row(i), {0, 0, 1}, exponentsMap[atomSym], pContractionMap[atomSym]};
            basisFunctionsList.push_back(AO_2px);
            basisFunctionsList.push_back(AO_2py);
            basisFunctionsList.push_back(AO_2pz);
        }
    }
    return basisFunctionsList;
}

/**
 * brief Print the molecule information
 * 
 */
void Molecule::printMoleculeInfo() {
    cout << "Number of atoms: " << nAtoms << endl;
    cout << "Charge: " << charge << endl;
    cout << "Atomic numbers: " << atomicNumbers.t();
    cout << "Atomic symbols: ";
    for (int i = 0; i < nAtoms; i++) {
        cout << atomicSymbols[i] << " ";
    }
    cout << endl;
    cout << "Coordinates:\n" << coordinates << endl;
    cout << "Number of basis functions: " << nBasisFunctions << endl;
    cout << "Number of electrons: " << nElectrons << endl;
    cout << "Number of alpha electrons: " << pAlpha << endl;
    cout << "Number of beta electrons: " << qBeta << endl;
}


/**
 * brief Calculate the analytical integration for the overlap of two primitive Gaussian shells for one dimension.
 * 
 * param alpha The exponent of the first primitive Gaussian shell
 * param beta The exponent of the second primitive Gaussian shell
 * param center_a The center of the first primitive Gaussian shell for one dimension
 * param center_b The center of the second primitive Gaussian shell for one dimension
 * param lA The angular momentum of the first primitive Gaussian shell
 * param lB The angular momentum of the second primitive Gaussian shell
 * 
 * return double The analytical overlap integral
 */
double overlapIntegral1D(double alpha, double beta, double center_a, double center_b, int lA, int lB) {
    // Calculate the exponential prefactor and the associated square root 
    double prefactor = exp(-alpha * beta * pow(center_a - center_b, 2) / (alpha + beta));
    prefactor *= sqrt(M_PI / (alpha + beta));

    // Calculate center_product
    double center_product = (alpha * center_a + beta * center_b) / (alpha + beta);

    // Double summation over the angular momentum combinations
    double sum = 0.0;
    for (int i = 0; i <= lA; i++) {
        for (int j = 0; j <= lB; j++) {
            // Only (i + j) even terms contribute
            if ((i + j) % 2 == 0) {
                sum += binomialCoef(lA, i) * binomialCoef(lB, j) * (doubleFactorial(i + j - 1) 
                        * pow(center_product - center_a, lA - i) * pow(center_product - center_b, lB - j)) 
                        / pow(2 * (alpha + beta), double(i + j) / 2);
            }
        }
    }
    double integral = prefactor * sum;
    return integral;
}

/**
 * brief Calculate the 3D overlap integral for two primitive Gaussian shells.
 * 
 * param centers_a The center of the first primitive Gaussian shell 
 * param centers_b The center of the second primitive Gaussian shell
 * param alpha The exponent of the first primitive Gaussian shell
 * param beta The exponent of the second primitive Gaussian shell
 * param lmn_a The angular momentum of the first primitive Gaussian shell
 * param lmn_b The angular momentum of the second primitive Gaussian shell
 * 
 * return double The overlap integral 
 */
double overlapIntegral3D(rowvec centers_a, rowvec centers_b, double alpha, double beta, ivec lmn_a, ivec lmn_b) {
    // Calculate the overlap integral for each dimension
    double integral = overlapIntegral1D(alpha, beta, centers_a(0), centers_b(0), lmn_a(0), lmn_b(0)) *
                      overlapIntegral1D(alpha, beta, centers_a(1), centers_b(1), lmn_a(1), lmn_b(1)) *
                      overlapIntegral1D(alpha, beta, centers_a(2), centers_b(2), lmn_a(2), lmn_b(2));

    return integral;
}

/**
 * 
 * @aram AO1 The first AO
 * param AO2 The second AO
 * 
 * return double The contracted overlap integral
 */
double calcContractedOverlap(AO AO1, AO AO2) {
    double contracted_overlap = 0.0;

    // Loop over all exponents
    for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
            double unnorm_overlap = overlapIntegral3D(AO1.center, AO2.center, 
                                                      AO1.exponents(k), AO2.exponents(l), 
                                                      AO1.lmn, AO2.lmn);
                                                      
            contracted_overlap += AO1.contraction_coeffs(k) * AO2.contraction_coeffs(l) *
                                  AO1.norm_constants(k) * AO2.norm_constants(l) * unnorm_overlap;
        }
    }
    return contracted_overlap;
}


/**
 * param molecule The molecule to calculate the overlap matrix for
 * 
 * return mat The contracted overlap matrix
 */
mat calcOverlapMatrix(Molecule molecule) {
    // Initialize the overlap matrix
    mat overlap_matrix = arma::zeros<mat>(molecule.nBasisFunctions, molecule.nBasisFunctions);

    // Loop over all basis function combinations
    for (int i = 0; i < molecule.nBasisFunctions; i++) {
        for (int j = 0; j < molecule.nBasisFunctions; j++) {
            overlap_matrix(i, j) = calcContractedOverlap(molecule.basisFunctionsList[i], molecule.basisFunctionsList[j]);
        }
    }

    return overlap_matrix;
}