#include <iostream>
#include "molecule.h"
#include "CNDO.h"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Incorrect file input! << endl";
        return 1;
    }

    const char* filename = argv[1];
    Molecule molecule(filename);

    // CNDO/2 SCF cycle
    mat S = calcOverlapMatrix(molecule);
    CNDO cndo(molecule, S);
    cndo.scfCycle();

    // Test X and Y matrices
    cout << "X Matrix: " << endl;
    cout << cndo.calcXMatrix() << endl;
    cout << "Y Matrix: " << endl;
    cout << cndo.calcYMatrix() << endl;

    // Test overlap gradient matrix
    cout << "Overlap Gradient Matrix: " << endl;
    cout << cndo.calcOverlapGradientMat() << endl;

    // Test gamma gradient matrix
    cout << "Gamma Gradient Matrix: " << endl;
    cout << cndo.calcGammaGradientMat() << endl;

    // Test nuclear repulsion energy gradient
    cout << "Nuclear Energy Gradient: " << endl;
    cout << cndo.calcNuclearEnergyGradMat() << endl;

    // Test total energy gradient
    cout << "Total Energy Gradient: " << endl;
    cout << cndo.calcTotalEnergyGradMat() << endl;

    return 0;
}