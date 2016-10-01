#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

/*
    P  - pressures,   atm;
    T  - temperature, celsius;
    sg - specific gravity (0.57 < sg < 1.68).
*/
double calcZfactor_DAK(double P, double T, double sg);

int main() {

    const double P  = 3250.0 * 6894.757293168 / 101325;
    const double T  = (213.0 - 32.0) * 5.0 / 9.0;
    const double sg = 0.666;
    const double z  = calcZfactor_DAK(P, T, sg);
    printf("Z = %f\n", z);

    return 0;
}

double calcZfactor_DAK(const double P, const double T, const double sg) {

    // Ppc - pseudocritical pressure, psia.
    // Tpc - pseudocritical temperature, K (degrees Rankine, 1(K) = 1*5/9 (°R)).
    // Suttons correlations, B.C. Craft and M.F. Hawkins.
    const double Ppc = 756.8 - 131.0 * sg - 3.60 * sg * sg;
    const double Tpc = (169.2 + 349.5 * sg - 74.0 * sg * sg) * 5.0 / 9.0;

    // Ppc - pseudo reduced pressure (1 (atm) = 1*101325/6894.757293168 (psia)).
    // Tpc - pseudo reduced temperature (1 (°C) = 1+273.15 (K)).
    // Dranchuk-Abbou Kassem: 0.2 < Ppc < 30, 1.0 < Tpc < 3.0.
    const double Ppr = P * 101325 / 6894.757293168 / Ppc;
    const double Tpr = (T + 273.15) / Tpc;
    printf("Ppr = %f, Tpr = %f\n", Ppr, Tpr);

    double invTpr[5];
    invTpr[0] = 1.0 / Tpr;
    for(uint8_t i = 1; i < 5; ++i)
        invTpr[i] = invTpr[i-1] * invTpr[0];

    double A[11];
    A[0] =  0.32650; A[1] = -1.07000; A[2] = -0.5339; A[3] = 0.01569;
    A[4] = -0.05165; A[5] =  0.54750; A[6] = -0.7361; A[7] = 0.18440;
    A[8] =  0.10560; A[9] =  0.61340; A[10] = 0.7210;

    double tmp;

    double C[4];
    C[0] = A[0] + A[1] * invTpr[0] + A[2] * invTpr[2] +
            A[3] * invTpr[3] + A[4] * invTpr[4];
    tmp = A[6] * invTpr[0] + A[7] * invTpr[1];
    C[1] = A[5] + tmp;
    C[2] = A[8] * tmp;

    double z = 1.0;
    uint16_t i = 0;
    uint16_t const maxLoop = 100;
    double convergence;

    // Newton's method
    for(i; i < maxLoop; ++i) {

        const double Rr  = 0.27*Ppr * invTpr[0] / z;
        const double Rr2 = Rr * Rr;
        const double Rr5 = Rr2 * Rr2 * Rr;

        tmp = A[10] * Rr2;

        C[3] = A[9] * (1.0 + tmp) * Rr2 * invTpr[2] * exp(-tmp);

        const double zn = 1.0 + C[0] * Rr + C[1] * Rr2 - C[2] * Rr5 + C[3];
        const double invZn = 1.0 / zn;
        const double dfdz = 1.0 + C[0] * Rr * invZn + 2.0*C[1] * Rr2 * invZn +
                             5.0*C[2] * Rr5 * invZn +
                             2.0*A[9] * Rr2 * invTpr[2] * invZn *
                              (1.0 + tmp - tmp * tmp * exp(-tmp));
        const double dz = (zn - z) / dfdz;

        z += dz;
        convergence = fabs(z - zn);
        if(convergence <= 1.0e-6)
            break;
    }

    printf("Loops = %" PRIu16 "\nConvergence = %e\n", i, convergence);
    if (i == maxLoop)
        printf("Warning: max loop = %" PRIu16 "\n", maxLoop);

    return z;

} const;